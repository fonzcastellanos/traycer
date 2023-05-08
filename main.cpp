#ifdef WIN32
#include <windows.h>
#endif

#if defined(WIN32) || defined(linux)
#include <GL/gl.h>
#include <GL/glut.h>
#elif defined(__APPLE__)
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#endif

#include <cassert>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <glm/glm.hpp>
#include <iostream>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "cli.hpp"
#include "main.hpp"
#include "stb_image_write.h"
#include "types.hpp"

#define IMG_W 640
#define IMG_H 480

#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100

const char *kWindowName = "traycer";

constexpr uint kImgArea = IMG_W * IMG_H;
constexpr float kAspectRatio = (float)IMG_W / IMG_H;

const float eps = 0.00000001;
const glm::vec3 backgroundColor(1.0, 1.0, 1.0);

static Config config;
static RenderMode render_mode = kRenderMode_Display;

// camera field of view
#define fov 60.0

// camera position
glm::vec3 camPos(0, 0, 0);

uchar buffer[IMG_H][IMG_W][3];

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
glm::vec3 ambient_light;

Light *extraLights[MAX_LIGHTS];

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g,
                        unsigned char b) {
  glColor3f(((float)r) / 255.0f, ((float)g) / 255.0f, ((float)b) / 255.0f);
  glVertex2i(x, y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g,
                     unsigned char b) {
  buffer[y][x][0] = r;
  buffer[y][x][1] = g;
  buffer[y][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g,
                unsigned char b) {
  plot_pixel_display(x, y, r, g, b);
  if (render_mode == kRenderMode_Jpeg) {
    plot_pixel_jpeg(x, y, r, g, b);
  }
}

void draw_scene() {
  for (unsigned int x = 0; x < IMG_W; x++) {
    glPointSize(2.0);
    glBegin(GL_POINTS);
    for (unsigned int y = 0; y < IMG_H; y++) {
      plot_pixel(x, y, buffer[y][x][0], buffer[y][x][1], buffer[y][x][2]);
    }
    glEnd();
    glFlush();
  }
  printf("Done!\n");
  fflush(stdout);
}

Status save_jpg() {
  printf("Saving JPEG file: %s\n", config.render_filepath);

  int rc = stbi_write_jpg(config.render_filepath, IMG_W, IMG_H,
                          kRgbChannel__Count, &buffer[0][0][0], 95);
  if (rc == 0) {
    std::fprintf(stderr, "Could not write data to JPEG file %s\n",
                 config.render_filepath);
    return kStatus_IoError;
  }

  return kStatus_Ok;
}

Status ParseField(FILE *f, const char *name, glm::vec3 *vals) {
  assert(f);
  assert(name);
  assert(vals);

  glm::vec3 &v_ = *vals;

  char str[100];

  int rc = std::fscanf(f, "%s", str);
  if (rc != 1) {
    std::fprintf(stderr,
                 "Failed to parse field name. Expected field name \"%s\".\n",
                 name);
    return kStatus_IoError;
  }

  if (std::strcmp(name, str) != 0) {
    std::fprintf(stderr,
                 "Received field name \"%s\". Expected field name \"%s\".\n",
                 str, name);
    // TODO: Choose a better-suited status.
    return kStatus_UnspecifiedError;
  };

  rc = std::fscanf(f, "%f %f %f", &v_[0], &v_[1], &v_[2]);
  if (rc != 3) {
    std::fprintf(stderr, "Failed to parse values of field \"%s\".\n", name);
    return kStatus_IoError;
  }

  std::printf("%s %f %f %f\n", name, v_[0], v_[1], v_[2]);

  return kStatus_Ok;
}

Status ParseField(std::FILE *f, const char *name, float *val) {
  assert(f);
  assert(name);
  assert(val);

  char str[100];

  int rc = std::fscanf(f, "%s", str);
  if (rc != 1) {
    std::fprintf(stderr,
                 "Failed to parse field name. Expected field name \"%s\".\n",
                 name);
    return kStatus_IoError;
  }

  if (std::strcmp(name, str) != 0) {
    std::fprintf(stderr,
                 "Received field name \"%s\". Expected field name \"%s\".\n",
                 str, name);
    // TODO: Choose a better-suited status.
    return kStatus_UnspecifiedError;
  };

  rc = std::fscanf(f, "%f", val);
  if (rc != 1) {
    std::fprintf(stderr, "Failed to parse value of field \"%s\".\n", name);
    return kStatus_IoError;
  }

  std::printf("%s: %f\n", name, *val);

  return kStatus_Ok;
}

Status ParseVertex(std::FILE *f, Vertex *v) {
  assert(f);
  assert(v);

  Status st = ParseField(f, "pos:", &v->position);
  if (st != kStatus_Ok) {
    std::fprintf(stderr, "Failed to parse field.\n");
    return st;
  }
  st = ParseField(f, "nor:", &v->normal);
  if (st != kStatus_Ok) {
    std::fprintf(stderr, "Failed to parse field.\n");
    return st;
  }
  st = ParseField(f, "dif:", &v->color_diffuse);
  if (st != kStatus_Ok) {
    std::fprintf(stderr, "Failed to parse field.\n");
    return st;
  }
  st = ParseField(f, "spe:", &v->color_specular);
  if (st != kStatus_Ok) {
    std::fprintf(stderr, "Failed to parse field.\n");
    return st;
  }

  st = ParseField(f, "shi:", &v->shininess);
  if (st != kStatus_Ok) {
    std::fprintf(stderr, "Failed to parse field.\n");
    return st;
  }

  return kStatus_Ok;
}

Status ParseTriangle(std::FILE *f, Triangle *t) {
  assert(f);
  assert(t);

  for (uint i = 0; i < 3; ++i) {
    Status st = ParseVertex(f, &t->v[i]);
    if (st != kStatus_Ok) {
      std::fprintf(stderr, "Failed to parse vertex.\n");
      return st;
    }
  }

  return kStatus_Ok;
}

Status ParseSphere(std::FILE *f, Sphere *s) {
  assert(f);
  assert(s);

  Status st = ParseField(f, "pos:", &s->position);
  if (st != kStatus_Ok) {
    std::fprintf(stderr, "Failed to parse field.\n");
    return st;
  }

  st = ParseField(f, "rad:", &s->radius);
  if (st != kStatus_Ok) {
    std::fprintf(stderr, "Failed to parse field.\n");
    return st;
  }

  st = ParseField(f, "dif:", &s->color_diffuse);
  if (st != kStatus_Ok) {
    std::fprintf(stderr, "Failed to parse field.\n");
    return st;
  }

  st = ParseField(f, "spe:", &s->color_specular);
  if (st != kStatus_Ok) {
    std::fprintf(stderr, "Failed to parse field.\n");
    return st;
  }

  st = ParseField(f, "shi:", &s->shininess);
  if (st != kStatus_Ok) {
    std::fprintf(stderr, "Failed to parse field.\n");
    return st;
  }

  return kStatus_Ok;
}

Status ParseLight(std::FILE *f, Light *l) {
  assert(f);
  assert(l);

  Status st = ParseField(f, "pos:", &l->position);
  if (st != kStatus_Ok) {
    std::fprintf(stderr, "Failed to parse field.\n");
    return st;
  }

  st = ParseField(f, "col:", &l->color);
  if (st != kStatus_Ok) {
    std::fprintf(stderr, "Failed to parse field.\n");
    return st;
  }

  return kStatus_Ok;
}

Status LoadScene(const char *filepath) {
  std::FILE *file = std::fopen(filepath, "r");
  if (!file) {
    std::fprintf(stderr, "Failed to open file %s.\n", filepath);
    return kStatus_IoError;
  }

  uint obj_count;
  char type[50];
  Triangle t;
  Sphere s;
  Light l;

  int rc = std::fscanf(file, "%u", &obj_count);
  if (rc != 1) {
    std::fprintf(stderr, "Failed to parse object count from file %s.\n",
                 filepath);
    return kStatus_IoError;
  }

  std::printf("Object count: %u\n", obj_count);

  Status st = ParseField(file, "amb:", &ambient_light);
  if (st != kStatus_Ok) {
    std::fprintf(stderr, "Failed to parse field of floats.\n");
    return kStatus_IoError;
  }

  for (uint i = 0; i < obj_count; ++i) {
    int rc = std::fscanf(file, "%s\n", type);
    if (rc != 1) {
      std::fprintf(stderr, "Failed to read object type.\n");
      return kStatus_IoError;
    }

    std::printf("Object type: %s\n", type);

    if (std::strcmp(type, "triangle") == 0) {
      st = ParseTriangle(file, &t);
      if (st != kStatus_Ok) {
        std::fprintf(stderr, "Failed to parse triangle.\n");
        return st;
      }

      if (num_triangles == MAX_TRIANGLES) {
        std::fprintf(stderr,
                     "Too many triangles. Increase MAX_TRIANGLES if more "
                     "triangles are desired.\n");
        // TODO: Choose a better-suited status.
        return kStatus_UnspecifiedError;
      }

      triangles[num_triangles] = t;
      ++num_triangles;
    } else if (std::strcmp(type, "sphere") == 0) {
      st = ParseSphere(file, &s);
      if (st != kStatus_Ok) {
        std::fprintf(stderr, "Failed to parse sphere.\n");
        return st;
      }

      if (num_spheres == MAX_SPHERES) {
        std::fprintf(stderr,
                     "Too many spheres. Increase MAX_SPHERES if more spheres "
                     "are desired.\n");
        // TODO: Choose a better-suited status.
        return kStatus_UnspecifiedError;
      }

      spheres[num_spheres] = s;
      ++num_spheres;
    } else if (std::strcmp(type, "light") == 0) {
      st = ParseLight(file, &l);
      if (st != kStatus_Ok) {
        std::fprintf(stderr, "Failed to parse light.\n");
        return st;
      }

      if (num_lights == MAX_LIGHTS) {
        std::fprintf(stderr,
                     "Too many lights. Increase MAX_LIGHTS if more lights are "
                     "desired.\n");
        // TODO: Choose a better-suited status.
        return kStatus_UnspecifiedError;
      }

      lights[num_lights] = l;
      ++num_lights;
    } else {
      std::fprintf(stderr, "Invalid object type \"%s\" in scene file %s\n",
                   type, filepath);

      // TODO: Choose a better-suited status.
      return kStatus_UnspecifiedError;
    }
  }

  return kStatus_Ok;
}

void display() {}

void init() {
  glMatrixMode(GL_PROJECTION);
  glOrtho(0, IMG_W, 0, IMG_H, 1, -1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0, 0, 0, 0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle() {
  // hack to make it only draw once
  static int once = 0;
  if (!once) {
    draw_scene();
    if (render_mode == kRenderMode_Jpeg) {
      Status status = save_jpg();
      if (status != kStatus_Ok) {
        std::fprintf(stderr, "Failed to save JPEG file %s.\n",
                     config.render_filepath);
#ifdef __APPLE__
        exit(EXIT_FAILURE);
#else
        glutLeaveMainLoop();
#endif
      }
    }
  }
  once = 1;
}

Ray **makeRays(SamplingMode mode) {
  // Generate image plane corners
  const float z = -1.0;
  float tanHalfFOV = glm::tan((fov / 2.0) * (M_PI / 180.0));
  glm::vec3 topRight(kAspectRatio * tanHalfFOV, tanHalfFOV, z);
  glm::vec3 topLeft(-kAspectRatio * tanHalfFOV, tanHalfFOV, z);
  glm::vec3 bottomRight(kAspectRatio * tanHalfFOV, -tanHalfFOV, z);
  glm::vec3 bottomLeft(-kAspectRatio * tanHalfFOV, -tanHalfFOV, z);

  // Allocate ray array
  Ray **rays = new Ray *[kImgArea];
  for (uint i = 0; i < kImgArea; ++i) {
    rays[i] = new Ray[config.rays_per_pixel];
  }

  // Calculate pixel dimensions
  float pWidth = (topRight.x - topLeft.x) / IMG_W;
  float pHeight = (topRight.y - bottomRight.y) / IMG_H;

  // Sample
  int rayInd = 0;
  if (mode == DEFAULT) {
    for (float y = bottomRight.y; y < topRight.y - pHeight / 2.0;
         y += pHeight) {
      for (float x = topLeft.x; x < topRight.x - pWidth / 2.0; x += pWidth) {
        rays[rayInd][0].position =
            glm::vec3(x + pWidth / 2.0, y + pHeight / 2.0, z);
        rays[rayInd][0].direction =
            glm::normalize(rays[rayInd][0].position - camPos);
        ++rayInd;
      }
    }
  } else if (mode == SUPER_JITTER) {
    std::srand(std::time(NULL));
    for (float y = bottomRight.y; y < topRight.y - pHeight / 2.0;
         y += pHeight) {
      for (float x = topLeft.x; x < topRight.x - pWidth / 2.0; x += pWidth) {
        for (int r = 0; r < config.rays_per_pixel; ++r) {
          float xOffset = pWidth * ((float)(std::rand()) / (float)(RAND_MAX));
          float yOffset = pHeight * ((float)(std::rand()) / (float)(RAND_MAX));
          rays[rayInd][r].position = glm::vec3(x + xOffset, y + yOffset, z);
          rays[rayInd][r].direction =
              glm::normalize(rays[rayInd][r].position - camPos);
        }
        ++rayInd;
      }
    }
  }
  return rays;
}

void getReflectedRay(Ray *ray, Intersection *in, Ray *reflectedRay) {
  glm::vec3 p = ray->position + in->t * ray->direction;
  glm::vec3 n;
  if (in->type == SPHERE) {
    int s = in->data.sphere.index;
    n = (p - spheres[s].position) / spheres[s].radius;
  } else if (in->type == TRIANGLE) {
    int tri = in->data.triangle.index;
    float alpha = in->data.triangle.alpha;
    float beta = in->data.triangle.beta;
    float gamma = in->data.triangle.gamma;

    n = triangles[tri].v[0].normal * alpha + triangles[tri].v[1].normal * beta +
        triangles[tri].v[2].normal * gamma;
    n = glm::normalize(n);
  }
  float vn = glm::clamp<float>(glm::dot(-ray->direction, n), 0.0, 1.0);
  reflectedRay->direction = glm::normalize(2.0f * vn * n + ray->direction);
  reflectedRay->position = p;
}

void sphereIntersect(Ray *ray, int s, Intersection *out) {
  if (ray == NULL || out == NULL) return;

  out->type = SPHERE;
  out->data.sphere.index = s;
  out->ray = ray;
  out->hit = false;

  // Calculate quadratic coefficients
  glm::vec3 centerToRay = ray->position - spheres[s].position;
  float b = 2.0 * glm::dot(ray->direction, centerToRay);
  float c = glm::pow(centerToRay.x, 2) + glm::pow(centerToRay.y, 2) +
            glm::pow(centerToRay.z, 2) - glm::pow(spheres[s].radius, 2);

  // Calculate discriminant
  float discriminant = glm::pow(b, 2) - 4.0 * c;

  // Calculate roots and determine t
  float t;
  if (discriminant < -eps) {  // misses sphere
    return;
  } else if (glm::abs(discriminant) < eps) {  // hits sphere at one point
    t = -b / 2.0;
  } else {  // hits sphere at two points
    float q = -0.5 * (b + glm::sign(b) * glm::sqrt(discriminant));
    float root1 = c / q;
    float root2 = q;
    t = glm::min(root1, root2);
  }

  if (t > eps) {  // hit because ahead of origin
    out->t = t;
    out->hit = true;
  }
}

void triangleIntersect(Ray *ray, int tri, Intersection *out) {
  if (ray == NULL || out == NULL) return;

  out->type = TRIANGLE;
  out->data.triangle.index = tri;
  out->ray = ray;
  out->hit = false;

  // Alias for vertices
  Vertex(&v)[3] = triangles[tri].v;

  // Calculate triangle normal
  glm::vec3 edge01 = v[1].position - v[0].position;
  glm::vec3 edge02 = v[2].position - v[0].position;
  glm::vec3 n = glm::cross(edge01, edge02);
  float twiceTriArea = glm::length(n);
  n = glm::normalize(n);

  // Check if plane is parallel to ray
  float nDotRay = glm::dot(n, ray->direction);
  if (glm::abs(nDotRay) < eps) return;

  // Check for intersection behind ray origin
  float d = -glm::dot(n, v[0].position);
  float t = -(glm::dot(n, ray->position) + d) / nDotRay;
  if (t < eps) return;

  glm::vec3 p = ray->position + t * ray->direction;

  glm::vec3 np;  // vector perpendicular to subtriangle

  // Inside-outside test
  glm::vec3 edge0p = p - v[0].position;
  np = glm::cross(edge01, edge0p);
  if (glm::dot(n, np) < -eps) return;

  float gamma = glm::length(np) / twiceTriArea;

  glm::vec3 edge12 = v[2].position - v[1].position;
  glm::vec3 edge1p = p - v[1].position;
  np = glm::cross(edge12, edge1p);
  if (glm::dot(n, np) < -eps) return;

  float alpha = glm::length(np) / twiceTriArea;

  glm::vec3 edge20 = v[0].position - v[2].position;
  glm::vec3 edge2p = p - v[2].position;
  np = glm::cross(edge20, edge2p);
  if (glm::dot(n, np) < -eps) return;

  float beta = glm::length(np) / twiceTriArea;

  out->data.triangle.alpha = alpha;
  out->data.triangle.beta = beta;
  out->data.triangle.gamma = gamma;
  out->t = t;
  out->hit = true;
}

void intersect(Ray *ray, Intersection *prev, Intersection *out) {
  if (ray == NULL || out == NULL) return;
  out->ray = ray;
  out->t = (float)FLT_MAX;
  out->hit = false;

  Intersection current;
  // sphere intersections
  for (int s = 0; s < num_spheres; ++s) {
    if (prev != NULL && prev->type == SPHERE && prev->data.sphere.index == s)
      continue;

    sphereIntersect(ray, s, &current);

    if (!current.hit || current.t > out->t + eps) continue;

    out->type = current.type;
    out->data.sphere.index = current.data.sphere.index;
    out->ray = current.ray;
    out->t = current.t;
    out->hit = current.hit;
  }
  // triangle intersections
  for (int tri = 0; tri < num_triangles; ++tri) {
    if (prev != NULL && prev->type == TRIANGLE &&
        prev->data.triangle.index == tri)
      continue;

    triangleIntersect(ray, tri, &current);

    if (!current.hit || current.t > out->t + eps) continue;

    out->type = current.type;
    out->data.triangle.index = current.data.triangle.index;
    out->data.triangle.alpha = current.data.triangle.alpha;
    out->data.triangle.beta = current.data.triangle.beta;
    out->data.triangle.gamma = current.data.triangle.gamma;
    out->ray = current.ray;
    out->t = current.t;
    out->hit = current.hit;
  }
}

float getPhongColor(float lightColor, float diffuse, float specular, float ln,
                    float rv, float shininess) {
  return lightColor * (diffuse * ln + specular * glm::pow(rv, shininess));
}

glm::vec3 shade(Intersection *surface, int bounces);

glm::vec3 traceRay(Ray *ray, Intersection *prev, int bounces) {
  Intersection closest;
  intersect(ray, prev, &closest);
  if (!closest.hit) {
    return backgroundColor;
  } else {
    return shade(&closest, bounces);
  }
}

glm::vec3 shade(Intersection *surface, int bounces) {
  int index;
  glm::vec3 p;  // point of intersection
  glm::vec3 n;  // normal;
  glm::vec3 color_diffuse;
  glm::vec3 color_specular;
  float shininess;

  // calculate surface-specific parameters
  if (surface->type == SPHERE) {
    index = surface->data.sphere.index;
    p = surface->ray->position +
        surface->t * (surface->ray->direction);  // point of intersection
    n = (p - spheres[index].position) / spheres[index].radius;  // normal
    color_diffuse = spheres[index].color_diffuse;
    color_specular = spheres[index].color_specular;
    shininess = spheres[index].shininess;
  } else if (surface->type == TRIANGLE) {
    index = surface->data.triangle.index;
    p = surface->ray->position + surface->t * (surface->ray->direction);

    float alpha = surface->data.triangle.alpha;
    float beta = surface->data.triangle.beta;
    float gamma = surface->data.triangle.gamma;

    n = triangles[index].v[0].normal * alpha +
        triangles[index].v[1].normal * beta +
        triangles[index].v[2].normal * gamma;
    n = glm::normalize(n);

    color_diffuse = triangles[index].v[0].color_diffuse * alpha +
                    triangles[index].v[1].color_diffuse * beta +
                    triangles[index].v[2].color_diffuse * gamma;

    color_specular = triangles[index].v[0].color_specular * alpha +
                     triangles[index].v[1].color_specular * beta +
                     triangles[index].v[2].color_specular * gamma;

    shininess = triangles[index].v[0].shininess * alpha +
                triangles[index].v[1].shininess * beta +
                triangles[index].v[2].shininess * gamma;
  }

  // launch shadow rays for local phong color
  Ray shadow;
  Intersection occluder;
  float toLight;
  glm::vec3 phongColor = ambient_light;
  float ln, rv;
  glm::vec3 reflection;
  for (int l = 0; l < num_lights; ++l) {
    // launch main shadow ray
    shadow.position = p;
    shadow.direction = glm::normalize(lights[l].position - p);

    intersect(&shadow, surface, &occluder);

    toLight = glm::length(lights[l].position - shadow.position);
    if (occluder.hit && occluder.t + eps < toLight) continue;

    ln = glm::clamp<float>(glm::dot(shadow.direction, n), 0.0, 1.0);
    reflection = glm::normalize(2.0f * ln * n - shadow.direction);
    rv = glm::clamp<float>(glm::dot(reflection, -surface->ray->direction), 0.0,
                           1.0);

    for (int c = 0; c < 3; ++c) {
      phongColor[c] += getPhongColor(
          lights[l].color[c] / (config.extra_lights_per_light + 1),
          color_diffuse[c], color_specular[c], ln, rv, shininess);
    }

    // launch extra shadow rays
    for (int e = 0; e < config.extra_lights_per_light; ++e) {
      shadow.position = p;
      shadow.direction = glm::normalize(extraLights[l][e].position - p);

      intersect(&shadow, surface, &occluder);

      toLight = glm::length(extraLights[l][e].position - shadow.position);
      if (occluder.hit && occluder.t + eps < toLight) continue;

      ln = glm::clamp<float>(glm::dot(shadow.direction, n), 0.0, 1.0);
      reflection = glm::normalize(2.0f * ln * n - shadow.direction);
      rv = glm::clamp<float>(glm::dot(reflection, -surface->ray->direction),
                             0.0, 1.0);

      for (int c = 0; c < 3; ++c) {
        phongColor[c] += getPhongColor(
            lights[l].color[c] / (config.extra_lights_per_light + 1),
            color_diffuse[c], color_specular[c], ln, rv, shininess);
      }
    }
  }

  // launch reflected rays
  if (bounces > 0) {
    Ray reflectedRay;
    getReflectedRay(surface->ray, surface, &reflectedRay);
    glm::vec3 reflectedColor = traceRay(&reflectedRay, surface, bounces - 1);
    for (int c = 0; c < 3; ++c) {
      phongColor[c] *= (1 - color_specular[c]);
      phongColor[c] += color_specular[c] * reflectedColor[c];
    }
  }

  return phongColor;
}

Status ParseConfig(uint argc, char *argv[], Config *c) {
  assert(argv);
  assert(c);

  cli::Opt opts[] = {
      {"rays-per-pixel", cli::kOptType_Int, &c->rays_per_pixel},
      {"reflection-bounces", cli::kOptType_Int, &c->reflection_bounces},
      {"extra-lights-per-light", cli::kOptType_Int, &c->extra_lights_per_light},
      {"render-file", cli::kOptType_String, c->render_filepath},
  };
  uint size = sizeof(opts) / sizeof(opts[0]);
  uint argi;
  cli::Status st = ParseOpts(argc, argv, opts, size, &argi);
  if (st != cli::kStatus_Ok) {
    std::fprintf(stderr, "Failed to parse options: %s\n",
                 cli::StatusMessage(st));
    return kStatus_UnspecifiedError;
  }
  if (argi >= argc) {
    std::fprintf(stderr, "Missing required scene file argument.\n");
    return kStatus_UnspecifiedError;
  }

  c->scene_filepath = argv[argi];

  return kStatus_Ok;
}

int main(int argc, char **argv) {
  glutInit(&argc, argv);

  Status st = ParseConfig(argc, argv, &config);
  if (st != kStatus_Ok) {
    std::fprintf(stderr, "Failed to parse config.\n");
    return EXIT_FAILURE;
  }

  if (config.render_filepath[0] != '\0') {
    render_mode = kRenderMode_Jpeg;
  }

  st = LoadScene(config.scene_filepath);
  if (st != kStatus_Ok) {
    std::fprintf(stderr, "Failed to load scene.\n");
    return EXIT_FAILURE;
  }

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0, 0);
  glutInitWindowSize(IMG_W, IMG_H);
  // TODO: Destroy window upon exiting program.
  glutCreateWindow(kWindowName);
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();

  if (config.extra_lights_per_light > 0) {
    std::srand(std::time(NULL));
    for (int l = 0; l < num_lights; ++l) {
      extraLights[l] = new Light[config.extra_lights_per_light];
      for (int e = 0; e < config.extra_lights_per_light; ++e) {
        float xOffset =
            (-1.0 + 2.0 * ((float)(std::rand()) / (float)(RAND_MAX)));
        float yOffset =
            (-1.0 + 2.0 * ((float)(std::rand()) / (float)(RAND_MAX)));
        // float zOffset = (-1.0
        // + 2.0*((float)(std::rand())/(float)(RAND_MAX)));
        glm::vec3 offset(xOffset, yOffset, 0.0);

        extraLights[l][e].position = lights[l].position + offset;
        extraLights[l][e].color = lights[l].color;
      }
    }
  }

  Ray **rays = NULL;
  if (config.rays_per_pixel == 1) {
    rays = makeRays(DEFAULT);
  } else if (config.rays_per_pixel > 1) {
    rays = makeRays(SUPER_JITTER);
  }

  for (uint p = 0; p < kImgArea; ++p) {
    uint x = p % IMG_W;
    uint y = p / IMG_W;
    glm::vec3 totalColor(0, 0, 0);
    for (int r = 0; r < config.rays_per_pixel; ++r) {
      totalColor += traceRay(&rays[p][r], NULL, config.reflection_bounces);
    }
    for (uint c = 0; c < 3; ++c) {
      buffer[y][x][c] =
          glm::clamp<float>(totalColor[c] / config.rays_per_pixel, 0, 1) * 255;
    }
  }
  delete[] rays;

  glutMainLoop();
}