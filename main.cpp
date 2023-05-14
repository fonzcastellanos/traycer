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
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <glm/glm.hpp>
#include <limits>
#include <random>
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "cli.hpp"
#include "main.hpp"
#include "scene.hpp"
#include "status.hpp"
#include "stb_image_write.h"
#include "types.hpp"

#define PI 3.14159

#define IMG_W 640
#define IMG_H 480

#define FOCAL_LEN 1
#define FOV 60

const char *kWindowName = "traycer";

constexpr uint kImgArea = IMG_W * IMG_H;

const float eps = 0.00000001;
const glm::vec3 backgroundColor(1, 1, 1);

static Config config;
static RenderTarget render_target = kRenderTarget_Window;

const glm::vec3 camera_position(0, 0, 0);

uchar buffer[IMG_H][IMG_W][3];

glm::vec3 ambient_light;

void PlotPixelDisplay(int x, int y, uchar r, uchar g, uchar b) {
  glColor3f(r * (1.0f / 255), g * (1.0f / 255), b * (1.0f / 255));
  glVertex2i(x, y);
}

void PlotPixelJpeg(int x, int y, uchar r, uchar g, uchar b) {
  buffer[y][x][0] = r;
  buffer[y][x][1] = g;
  buffer[y][x][2] = b;
}

void PlotPixel(int x, int y, uchar r, uchar g, uchar b) {
  PlotPixelDisplay(x, y, r, g, b);
  if (render_target == kRenderTarget_Jpeg) {
    PlotPixelJpeg(x, y, r, g, b);
  }
}

void DrawScene() {
  for (uint x = 0; x < IMG_W; ++x) {
    glPointSize(2);
    glBegin(GL_POINTS);
    for (uint y = 0; y < IMG_H; ++y) {
      PlotPixel(x, y, buffer[y][x][0], buffer[y][x][1], buffer[y][x][2]);
    }
    glEnd();
    glFlush();
  }
  std::printf("Done!\n");
  fflush(stdout);
}

Status SaveJpeg() {
  std::printf("Saving JPEG file %s.\n", config.render_filepath);

  stbi_flip_vertically_on_write(1);
  int rc = stbi_write_jpg(config.render_filepath, IMG_W, IMG_H,
                          kRgbChannel__Count, &buffer[0][0][0], 95);
  if (rc == 0) {
    std::fprintf(stderr, "Could not write data to JPEG file %s.\n",
                 config.render_filepath);
    return kStatus_IoError;
  }

  return kStatus_Ok;
}

void Display() {}

void Init() {
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
    DrawScene();
    if (render_target == kRenderTarget_Jpeg) {
      Status status = SaveJpeg();
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

Ray **MakeRays(SamplingMode mode, uint w, uint h, float fov, float focal_len) {
  // Calculate projection plane corners.
  float proj_plane_z = -focal_len;
  float aspect = (float)w / h;
  float tan_half_fov = glm::tan((fov / 2) * ((float)PI / 180));
  glm::vec3 tr(aspect * tan_half_fov, tan_half_fov, proj_plane_z);
  glm::vec3 tl(-aspect * tan_half_fov, tan_half_fov, proj_plane_z);
  glm::vec3 br(aspect * tan_half_fov, -tan_half_fov, proj_plane_z);
  glm::vec3 bl(-aspect * tan_half_fov, -tan_half_fov, proj_plane_z);

  constexpr float kTolerance = 0.00001;

  uint pixel_count = w * h;

  Ray **rays = new Ray *[pixel_count];
  for (uint i = 0; i < pixel_count; ++i) {
    rays[i] = new Ray[config.rays_per_pixel];
  }

  float pix_w = (tr.x - tl.x) / w;
  float pix_h = (tr.y - br.y) / h;

  switch (mode) {
    case kSamplingMode_Default: {
      uint ray_idx = 0;
      for (float y = br.y + pix_h * (1.0f / 2); y < tr.y; y += pix_h) {
        for (float x = tl.x + pix_w * (1.0f / 2); x < tr.x; x += pix_w) {
          rays[ray_idx][0].position = glm::vec3(x, y, proj_plane_z);
          rays[ray_idx][0].direction =
              glm::normalize(rays[ray_idx][0].position - camera_position);
          ++ray_idx;
        }
      }

      break;
    }
    case kSamplingMode_Jitter: {
      std::default_random_engine re(std::random_device{}());
      std::uniform_real_distribution<float> distrib(0, 1);

      uint ray_idx = 0;
      for (float y = br.y; y < tr.y - kTolerance; y += pix_h) {
        for (float x = tl.x; x < tr.x - kTolerance; x += pix_w) {
          for (int i = 0; i < config.rays_per_pixel; ++i) {
            float x_offset = pix_w * distrib(re);
            float y_offset = pix_h * distrib(re);
            rays[ray_idx][i].position =
                glm::vec3(x + x_offset, y + y_offset, proj_plane_z);
            rays[ray_idx][i].direction =
                glm::normalize(rays[ray_idx][i].position - camera_position);
          }
          ++ray_idx;
        }
      }

      break;
    }
  }

  return rays;
}

void GetReflectedRay(Ray *ray, Sphere *spheres, Triangle *triangles,
                     Intersection *in, Ray *reflected_ray) {
  assert(ray);
  assert(spheres);
  assert(triangles);
  assert(in);
  assert(reflected_ray);

  glm::vec3 p = ray->position + in->t * ray->direction;
  glm::vec3 n;
  if (in->type == kGeometryType_Sphere) {
    int s = in->sphere_data.index;
    n = (p - spheres[s].position) / spheres[s].radius;
  } else if (in->type == kGeometryType_Triangle) {
    int tri = in->triangle_data.index;
    float alpha = in->triangle_data.alpha;
    float beta = in->triangle_data.beta;
    float gamma = in->triangle_data.gamma;

    n = triangles[tri].v[0].normal * alpha + triangles[tri].v[1].normal * beta +
        triangles[tri].v[2].normal * gamma;
    n = glm::normalize(n);
  }
  float vn = glm::clamp<float>(glm::dot(-ray->direction, n), 0.0, 1.0);
  reflected_ray->direction = glm::normalize(2.0f * vn * n + ray->direction);
  reflected_ray->position = p;
}

void IntersectSphere(Ray *ray, Sphere *spheres, uint sphere_idx,
                     Intersection *out) {
  assert(ray);
  assert(spheres);
  assert(out);

  out->type = kGeometryType_Sphere;
  out->sphere_data.index = sphere_idx;
  out->ray = ray;
  out->hit = false;

  // Calculate quadratic coefficients
  glm::vec3 center_to_ray = ray->position - spheres[sphere_idx].position;
  float b = 2 * glm::dot(ray->direction, center_to_ray);
  float c = glm::pow(center_to_ray.x, 2) + glm::pow(center_to_ray.y, 2) +
            glm::pow(center_to_ray.z, 2) -
            glm::pow(spheres[sphere_idx].radius, 2);

  // Calculate discriminant
  float discriminant = glm::pow(b, 2) - 4 * c;

  // Calculate roots and determine t
  float t;
  if (discriminant < -eps) {  // misses sphere
    return;
  } else if (glm::abs(discriminant) < eps) {  // hits sphere at one point
    t = -b * (1.0f / 2);
  } else {  // hits sphere at two points
    float q = (-1.0f / 2) * (b + glm::sign(b) * glm::sqrt(discriminant));
    float root1 = c / q;
    float root2 = q;
    t = glm::min(root1, root2);
  }

  if (t > eps) {  // hit because ahead of origin
    out->t = t;
    out->hit = true;
  }
}

void IntersectTriangle(Ray *ray, Triangle *triangles, int triangle_idx,
                       Intersection *out) {
  assert(ray);
  assert(out);

  out->type = kGeometryType_Triangle;
  out->triangle_data.index = triangle_idx;
  out->ray = ray;
  out->hit = false;

  // Alias for vertices
  Vertex(&v)[3] = triangles[triangle_idx].v;

  // Calculate triangle normal
  glm::vec3 edge01 = v[1].position - v[0].position;
  glm::vec3 edge02 = v[2].position - v[0].position;
  glm::vec3 n = glm::cross(edge01, edge02);
  float twiceTriArea = glm::length(n);
  n = glm::normalize(n);

  // Check if plane is parallel to ray
  float nDotRay = glm::dot(n, ray->direction);
  if (glm::abs(nDotRay) < eps) {
    return;
  }

  // Check for intersection behind ray origin
  float d = -glm::dot(n, v[0].position);
  float t = -(glm::dot(n, ray->position) + d) / nDotRay;
  if (t < eps) {
    return;
  }

  glm::vec3 p = ray->position + t * ray->direction;

  glm::vec3 np;  // vector perpendicular to subtriangle

  // Inside-outside test
  glm::vec3 edge0p = p - v[0].position;
  np = glm::cross(edge01, edge0p);
  if (glm::dot(n, np) < -eps) {
    return;
  }

  float gamma = glm::length(np) / twiceTriArea;

  glm::vec3 edge12 = v[2].position - v[1].position;
  glm::vec3 edge1p = p - v[1].position;
  np = glm::cross(edge12, edge1p);
  if (glm::dot(n, np) < -eps) {
    return;
  }

  float alpha = glm::length(np) / twiceTriArea;

  glm::vec3 edge20 = v[0].position - v[2].position;
  glm::vec3 edge2p = p - v[2].position;
  np = glm::cross(edge20, edge2p);
  if (glm::dot(n, np) < -eps) {
    return;
  }

  float beta = glm::length(np) / twiceTriArea;

  out->triangle_data.alpha = alpha;
  out->triangle_data.beta = beta;
  out->triangle_data.gamma = gamma;
  out->t = t;
  out->hit = true;
}

void Intersect(Ray *ray, Sphere *spheres, uint sphere_count,
               Triangle *triangles, uint triangle_count, Intersection *prev,
               Intersection *out) {
  assert(ray);
  assert(spheres);
  assert(triangles);
  assert(prev);
  assert(out);

  out->ray = ray;
  out->t = std::numeric_limits<float>::max();
  out->hit = false;

  Intersection current;
  // sphere intersections
  for (uint i = 0; i < sphere_count; ++i) {
    if (prev != NULL && prev->type == kGeometryType_Sphere &&
        prev->sphere_data.index == i) {
      continue;
    }

    IntersectSphere(ray, spheres, i, &current);

    if (!current.hit || current.t > out->t + eps) {
      continue;
    }

    out->type = current.type;
    out->sphere_data.index = current.sphere_data.index;
    out->ray = current.ray;
    out->t = current.t;
    out->hit = current.hit;
  }
  // triangle intersections
  for (uint i = 0; i < triangle_count; ++i) {
    if (prev != NULL && prev->type == kGeometryType_Triangle &&
        prev->triangle_data.index == i) {
      continue;
    }

    IntersectTriangle(ray, triangles, i, &current);

    if (!current.hit || current.t > out->t + eps) {
      continue;
    }

    out->type = current.type;
    out->triangle_data.index = current.triangle_data.index;
    out->triangle_data.alpha = current.triangle_data.alpha;
    out->triangle_data.beta = current.triangle_data.beta;
    out->triangle_data.gamma = current.triangle_data.gamma;
    out->ray = current.ray;
    out->t = current.t;
    out->hit = current.hit;
  }
}

float GetPhongColor(float lightColor, float diffuse, float specular, float ln,
                    float rv, float shininess) {
  return lightColor * (diffuse * ln + specular * glm::pow(rv, shininess));
}

glm::vec3 Shade(Intersection *surface, Sphere *spheres, uint sphere_count,
                Triangle *triangles, uint triangle_count, Light *lights,
                uint light_count, int bounces, const Lights *extra_lights);

glm::vec3 TraceRay(Ray *ray, Sphere *spheres, uint sphere_count,
                   Triangle *triangles, uint triangle_count, Light *lights,
                   uint light_count, Intersection *prev, int bounces,
                   const Lights *extra_lights) {
  Intersection closest;
  Intersect(ray, spheres, sphere_count, triangles, triangle_count, prev,
            &closest);
  if (!closest.hit) {
    return backgroundColor;
  } else {
    return Shade(&closest, spheres, sphere_count, triangles, triangle_count,
                 lights, light_count, bounces, extra_lights);
  }
}

glm::vec3 Shade(Intersection *surface, Sphere *spheres, uint sphere_count,
                Triangle *triangles, uint triangle_count, Light *lights,
                uint light_count, int bounces, const Lights *extra_lights) {
  assert(surface);
  assert(spheres);
  assert(triangles);
  assert(lights);

  int index;
  glm::vec3 p;  // point of intersection
  glm::vec3 n;  // normal;
  glm::vec3 color_diffuse;
  glm::vec3 color_specular;
  float shininess;

  // calculate surface-specific parameters
  if (surface->type == kGeometryType_Sphere) {
    index = surface->sphere_data.index;
    p = surface->ray->position +
        surface->t * (surface->ray->direction);  // point of intersection
    n = (p - spheres[index].position) / spheres[index].radius;  // normal
    color_diffuse = spheres[index].color_diffuse;
    color_specular = spheres[index].color_specular;
    shininess = spheres[index].shininess;
  } else if (surface->type == kGeometryType_Triangle) {
    index = surface->triangle_data.index;
    p = surface->ray->position + surface->t * (surface->ray->direction);

    float alpha = surface->triangle_data.alpha;
    float beta = surface->triangle_data.beta;
    float gamma = surface->triangle_data.gamma;

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
  for (uint l = 0; l < light_count; ++l) {
    // launch main shadow ray
    shadow.position = p;
    shadow.direction = glm::normalize(lights[l].position - p);

    Intersect(&shadow, spheres, sphere_count, triangles, triangle_count,
              surface, &occluder);

    toLight = glm::length(lights[l].position - shadow.position);
    if (occluder.hit && occluder.t + eps < toLight) {
      continue;
    }

    ln = glm::clamp<float>(glm::dot(shadow.direction, n), 0.0, 1.0);
    reflection = glm::normalize(2.0f * ln * n - shadow.direction);
    rv = glm::clamp<float>(glm::dot(reflection, -surface->ray->direction), 0.0,
                           1.0);

    for (int c = 0; c < 3; ++c) {
      phongColor[c] += GetPhongColor(
          lights[l].color[c] / (config.extra_lights_per_light + 1),
          color_diffuse[c], color_specular[c], ln, rv, shininess);
    }

    // launch extra shadow rays
    for (uint e = 0; e < config.extra_lights_per_light; ++e) {
      shadow.position = p;
      uint k = l * config.extra_lights_per_light + e;
      shadow.direction = glm::normalize(extra_lights->positions[k] - p);

      Intersect(&shadow, spheres, sphere_count, triangles, triangle_count,
                surface, &occluder);

      toLight = glm::length(extra_lights->positions[k] - shadow.position);
      if (occluder.hit && occluder.t + eps < toLight) {
        continue;
      }

      ln = glm::clamp<float>(glm::dot(shadow.direction, n), 0.0, 1.0);
      reflection = glm::normalize(2.0f * ln * n - shadow.direction);
      rv = glm::clamp<float>(glm::dot(reflection, -surface->ray->direction),
                             0.0, 1.0);

      for (int c = 0; c < 3; ++c) {
        phongColor[c] += GetPhongColor(
            lights[l].color[c] / (config.extra_lights_per_light + 1),
            color_diffuse[c], color_specular[c], ln, rv, shininess);
      }
    }
  }

  // launch reflected rays
  if (bounces > 0) {
    Ray reflectedRay;
    GetReflectedRay(surface->ray, spheres, triangles, surface, &reflectedRay);
    glm::vec3 reflectedColor = TraceRay(
        &reflectedRay, spheres, sphere_count, triangles, triangle_count, lights,
        light_count, surface, bounces - 1, extra_lights);
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
      {"extra-lights-per-light", cli::kOptType_Uint,
       &c->extra_lights_per_light},
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
    render_target = kRenderTarget_Jpeg;
  }

  std::vector<Sphere> spheres;
  std::vector<Triangle> triangles;
  std::vector<Light> lights;

  st = LoadScene(config.scene_filepath, &ambient_light, &triangles, &spheres,
                 &lights);
  if (st != kStatus_Ok) {
    std::fprintf(stderr, "Failed to load scene.\n");
    return EXIT_FAILURE;
  }

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0, 0);
  glutInitWindowSize(IMG_W, IMG_H);
  // TODO: Destroy window upon exiting program.
  glutCreateWindow(kWindowName);
  glutDisplayFunc(Display);
  glutIdleFunc(idle);
  Init();

  Lights extra_lights;

  if (config.extra_lights_per_light > 0) {
    std::default_random_engine re(std::random_device{}());
    std::uniform_real_distribution<float> distrib(-1, 1);

    uint count = lights.size() * config.extra_lights_per_light;
    extra_lights.positions.resize(count);
    extra_lights.colors.resize(count);

    for (uint i = 0; i < lights.size(); ++i) {
      for (uint j = 0; j < config.extra_lights_per_light; ++j) {
        float x = distrib(re);
        float y = distrib(re);
        glm::vec3 offset(x, y, 0);

        uint k = i * config.extra_lights_per_light + j;

        extra_lights.positions[k] = lights[i].position + offset;
        extra_lights.colors[k] = lights[i].color;
      }
    }
  }

  Ray **rays = 0;
  if (config.rays_per_pixel == 1) {
    rays = MakeRays(kSamplingMode_Default, IMG_W, IMG_H, FOV, FOCAL_LEN);
  } else if (config.rays_per_pixel > 1) {
    rays = MakeRays(kSamplingMode_Jitter, IMG_W, IMG_H, FOV, FOCAL_LEN);
  }

  for (uint i = 0; i < kImgArea; ++i) {
    uint x = i % IMG_W;
    uint y = i / IMG_W;
    glm::vec3 total_color(0, 0, 0);
    for (int j = 0; j < config.rays_per_pixel; ++j) {
      total_color += TraceRay(&rays[i][j], spheres.data(), spheres.size(),
                              triangles.data(), triangles.size(), lights.data(),
                              lights.size(), NULL, config.reflection_bounces,
                              &extra_lights);
    }
    for (int j = 0; j < kRgbChannel__Count; ++j) {
      buffer[y][x][j] =
          glm::clamp<float>(total_color[j] / config.rays_per_pixel, 0, 1) * 255;
    }
  }
  delete[] rays;

  glutMainLoop();
}