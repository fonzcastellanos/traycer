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

constexpr float kTolerance = 0.000001;
constexpr float kShadowRayBias = kTolerance;
constexpr float kReflectionRayBias = kTolerance;

const glm::vec3 kBackgroundColor(1, 1, 1);
const glm::vec3 kCameraPosition(0, 0, 0);

static Config config;
static RenderTarget render_target = kRenderTarget_Window;

uchar buffer[IMG_H][IMG_W][3];

glm::vec3 ambient_light;

static void PlotPixelDisplay(int x, int y, uchar r, uchar g, uchar b) {
  glColor3f(r * (1.0f / 255), g * (1.0f / 255), b * (1.0f / 255));
  glVertex2i(x, y);
}

static void PlotPixelJpeg(int x, int y, uchar r, uchar g, uchar b) {
  buffer[y][x][0] = r;
  buffer[y][x][1] = g;
  buffer[y][x][2] = b;
}

static void PlotPixel(int x, int y, uchar r, uchar g, uchar b) {
  PlotPixelDisplay(x, y, r, g, b);
  if (render_target == kRenderTarget_Jpeg) {
    PlotPixelJpeg(x, y, r, g, b);
  }
}

static void DrawScene() {
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

static Status SaveJpeg() {
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

static void Display() {}

static void Init() {
  glMatrixMode(GL_PROJECTION);
  glOrtho(0, IMG_W, 0, IMG_H, 1, -1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0, 0, 0, 0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void Idle() {
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

static void ProjectionPlaneCorners(uint w, uint h, float fov, float z,
                                   glm::vec3 *tr, glm::vec3 *tl, glm::vec3 *br,
                                   glm::vec3 *bl) {
  assert(tr);
  assert(tl);
  assert(br);
  assert(bl);

  assert(fov > kTolerance);

  float aspect = (float)w / h;
  float tan_half_fov = glm::tan(fov * 0.5f * (float)PI * (1.0f / 180));

  *tr = glm::vec3(aspect * tan_half_fov, tan_half_fov, z);
  *tl = glm::vec3(-aspect * tan_half_fov, tan_half_fov, z);
  *br = glm::vec3(aspect * tan_half_fov, -tan_half_fov, z);
  *bl = glm::vec3(-aspect * tan_half_fov, -tan_half_fov, z);
}

static void MakeDefaultRays(const glm::vec3 *camera_position, uint w, uint h,
                            float fov, float focal_len,
                            std::vector<Ray> *rays) {
  assert(camera_position);
  assert(fov + kTolerance > 0);
  assert(focal_len + kTolerance > 0);
  assert(rays);

  float projection_plane_z = -focal_len;
  glm::vec3 tr;
  glm::vec3 tl;
  glm::vec3 br;
  glm::vec3 bl;
  ProjectionPlaneCorners(w, h, fov, projection_plane_z, &tr, &tl, &br, &bl);

  uint pixel_count = w * h;

  auto &rays_ = *rays;

  rays_.resize(pixel_count);

  float pix_w = (tr.x - tl.x) / w;
  float pix_h = (tr.y - br.y) / h;

  uint i = 0;
  for (float y = br.y + pix_h * 0.5f; y < tr.y; y += pix_h) {
    for (float x = tl.x + pix_w * 0.5f; x < tr.x; x += pix_w) {
      rays_[i].position = glm::vec3(x, y, projection_plane_z);
      rays_[i].direction = glm::normalize(rays_[i].position - *camera_position);
      ++i;
    }
  }
}

static void MakeJitteredRays(const glm::vec3 *camera_position, uint w, uint h,
                             float fov, float focal_len, uint rays_per_pixel,
                             std::vector<Ray> *rays) {
  assert(camera_position);
  assert(fov + kTolerance > 0);
  assert(focal_len + kTolerance > 0);
  assert(rays_per_pixel > 0);
  assert(rays);

  float projection_plane_z = -focal_len;
  glm::vec3 tr;
  glm::vec3 tl;
  glm::vec3 br;
  glm::vec3 bl;
  ProjectionPlaneCorners(w, h, fov, projection_plane_z, &tr, &tl, &br, &bl);

  uint pixel_count = w * h;

  auto &rays_ = *rays;

  rays_.resize(pixel_count * rays_per_pixel);

  float pix_w = (tr.x - tl.x) / w;
  float pix_h = (tr.y - br.y) / h;

  std::default_random_engine re(std::random_device{}());
  std::uniform_real_distribution<float> distrib(0, 1);

  uint i = 0;
  for (float y = br.y; y < tr.y - kTolerance; y += pix_h) {
    for (float x = tl.x; x < tr.x - kTolerance; x += pix_w) {
      for (uint j = 0; j < rays_per_pixel; ++j) {
        float x_offset = pix_w * distrib(re);
        float y_offset = pix_h * distrib(re);

        uint k = i * rays_per_pixel + j;

        rays_[k].position =
            glm::vec3(x + x_offset, y + y_offset, projection_plane_z);
        rays_[k].direction =
            glm::normalize(rays_[k].position - *camera_position);
      }
      ++i;
    }
  }
}

static float IntersectSphere(Ray *ray, const Sphere *sph) {
  assert(ray);
  assert(sph);

  // Calculate quadratic coefficients.
  glm::vec3 sph_to_ray = ray->position - sph->position;
  float b = 2 * glm::dot(ray->direction, sph_to_ray);
  float c = sph_to_ray.x * sph_to_ray.x + sph_to_ray.y * sph_to_ray.y +
            sph_to_ray.z * sph_to_ray.z - sph->radius * sph->radius;

  float discriminant = b * b - 4 * c;

  float t;
  if (discriminant < -kTolerance) {
    // Misses sphere.
    t = 0;
  } else if (discriminant < kTolerance) {
    // Hits sphere at one point.
    t = -b * 0.5f;
  } else {
    // Hits sphere at two points.

    // Implemented this way to improve stability of calculation. Addresses the
    // issue of when b and the root of the discriminant don't have the same sign
    // but the values are very close to each other.
    float q = -0.5f * (b + glm::sign(b) * glm::sqrt(discriminant));
    float root1 = c / q;
    float root2 = q;
    t = glm::min(root1, root2);
  }

  return t;
}

static float IntersectTriangle(Ray *ray, const Triangle *tri, float *alpha,
                               float *beta, float *gamma) {
  assert(ray);
  assert(tri);
  assert(alpha);
  assert(beta);
  assert(gamma);

  const auto &verts = tri->vertices;

  float n_dot_ray = glm::dot(tri->normal, ray->direction);

  int plane_parallel_to_ray = glm::abs(n_dot_ray) < kTolerance;
  if (plane_parallel_to_ray) {
    return 0;
  }

  float d = -glm::dot(tri->normal, verts[0].position);
  float t = -(glm::dot(tri->normal, ray->position) + d) / n_dot_ray;

  int intersection_behind_ray_origin = t < kTolerance;
  if (intersection_behind_ray_origin) {
    return 0;
  }

  glm::vec3 p = ray->position + t * ray->direction;

  // Perpendicular to subtriangle.
  glm::vec3 np;

  glm::vec3 e0p = p - verts[0].position;
  glm::vec3 e01 = verts[1].position - verts[0].position;
  np = glm::cross(e01, e0p);
  // Inside-outside test.
  if (glm::dot(tri->normal, np) < -kTolerance) {
    return 0;
  }

  float twice_area = 2 * tri->area;

  float gamma_ = glm::length(np) / twice_area;

  glm::vec3 e12 = verts[2].position - verts[1].position;
  glm::vec3 e1p = p - verts[1].position;
  np = glm::cross(e12, e1p);
  if (glm::dot(tri->normal, np) < -kTolerance) {
    return 0;
  }

  float alpha_ = glm::length(np) / twice_area;

  glm::vec3 e20 = verts[0].position - verts[2].position;
  glm::vec3 e2p = p - verts[2].position;
  np = glm::cross(e20, e2p);
  if (glm::dot(tri->normal, np) < -kTolerance) {
    return 0;
  }

  *alpha = alpha_;
  *gamma = gamma_;
  *beta = 1 - alpha_ - gamma_;

  return t;
}

static int Intersect(Ray *ray, const Sphere *spheres, uint sphere_count,
                     const Triangle *triangles, uint triangle_count,
                     Intersection *intersection) {
  assert(ray);
  assert(spheres);
  assert(triangles);
  assert(intersection);

  intersection->ray = ray;
  intersection->t = std::numeric_limits<float>::max() - kTolerance;

  int found_intersection = 0;

  // Spheres
  for (uint i = 0; i < sphere_count; ++i) {
    float t = IntersectSphere(ray, &spheres[i]);

    int hit = t > kTolerance;
    int further = t > intersection->t + kTolerance;
    if (!hit || further) {
      continue;
    }

    intersection->type = kGeometryType_Sphere;
    intersection->sphere.index = i;
    intersection->t = t;
    found_intersection = 1;
  }

  // Triangles
  for (uint i = 0; i < triangle_count; ++i) {
    float alpha;
    float beta;
    float gamma;
    float t = IntersectTriangle(ray, &triangles[i], &alpha, &beta, &gamma);

    int hit = t > kTolerance;
    int further = t > intersection->t + kTolerance;
    if (!hit || further) {
      continue;
    }

    intersection->type = kGeometryType_Triangle;
    intersection->triangle.index = i;
    intersection->triangle.alpha = alpha;
    intersection->triangle.beta = beta;
    intersection->triangle.gamma = gamma;
    intersection->t = t;
    found_intersection = 1;
  }

  return found_intersection;
}

static glm::vec3 PhongColor(glm::vec3 light_color, glm::vec3 diffuse_color,
                            glm::vec3 specular_color, float l_dot_n,
                            float r_dot_v, float shininess) {
  glm::vec3 result =
      light_color * ((diffuse_color * l_dot_n) +
                     (specular_color * glm::pow(r_dot_v, shininess)));
  return result;
}

static glm::vec3 Shade(Intersection *surface, const Scene *scene, int bounces,
                       const Lights *extra_lights);

static glm::vec3 TraceRay(Ray *ray, const Scene *scene, int bounces,
                          const Lights *extra_lights) {
  assert(ray);
  assert(scene);

  Intersection intersection;
  int hit = Intersect(ray, scene->spheres, scene->sphere_count,
                      scene->triangles, scene->triangle_count, &intersection);
  if (!hit) {
    return kBackgroundColor;
  } else {
    return Shade(&intersection, scene, bounces, extra_lights);
  }
}

static glm::vec3 Shade(Intersection *surface, const Scene *scene, int bounces,
                       const Lights *extra_lights) {
  assert(surface);
  assert(scene);
  assert(extra_lights);

  const Triangle *triangles = scene->triangles;
  uint triangle_count = scene->triangle_count;
  const Sphere *spheres = scene->spheres;
  uint sphere_count = scene->sphere_count;
  const Light *lights = scene->lights;
  uint light_count = scene->light_count;

  glm::vec3 intersection_position =
      surface->ray->position + surface->t * surface->ray->direction;
  glm::vec3 normal;
  glm::vec3 color_diffuse;
  glm::vec3 color_specular;
  float shininess;

  // Calculate surface-specific parameters.
  if (surface->type == kGeometryType_Sphere) {
    const auto &sph = spheres[surface->sphere.index];

    normal = (intersection_position - sph.position) / sph.radius;

    color_diffuse = sph.color_diffuse;
    color_specular = sph.color_specular;
    shininess = sph.shininess;
  } else if (surface->type == kGeometryType_Triangle) {
    const auto &tri = triangles[surface->triangle.index];

    float alpha = surface->triangle.alpha;
    float beta = surface->triangle.beta;
    float gamma = surface->triangle.gamma;

    normal = tri.vertices[0].normal * alpha + tri.vertices[1].normal * beta +
             tri.vertices[2].normal * gamma;
    normal = glm::normalize(normal);

    color_diffuse = tri.vertices[0].color_diffuse * alpha +
                    tri.vertices[1].color_diffuse * beta +
                    tri.vertices[2].color_diffuse * gamma;

    color_specular = tri.vertices[0].color_specular * alpha +
                     tri.vertices[1].color_specular * beta +
                     tri.vertices[2].color_specular * gamma;

    shininess = tri.vertices[0].shininess * alpha +
                tri.vertices[1].shininess * beta +
                tri.vertices[2].shininess * gamma;
  }

  // Launch shadow rays for local phong color.
  Intersection occluder;
  glm::vec3 phong_color = ambient_light;
  float l_dot_n;
  float r_dot_v;
  glm::vec3 reflection;
  for (uint l = 0; l < light_count; ++l) {
    // Launch main shadow ray.
    Ray main_shadow_ray;
    main_shadow_ray.position = intersection_position + kShadowRayBias * normal;

    glm::vec3 shadow_to_light = lights[l].position - main_shadow_ray.position;
    float shadow_to_light_dist = glm::length(shadow_to_light);

    main_shadow_ray.direction = shadow_to_light / shadow_to_light_dist;

    int hit = Intersect(&main_shadow_ray, spheres, sphere_count, triangles,
                        triangle_count, &occluder);

    int occluded = hit && occluder.t + kTolerance < shadow_to_light_dist;
    if (occluded) {
      continue;
    }

    l_dot_n =
        glm::clamp<float>(glm::dot(main_shadow_ray.direction, normal), 0, 1);
    reflection =
        glm::normalize(2 * l_dot_n * normal - main_shadow_ray.direction);
    r_dot_v =
        glm::clamp<float>(glm::dot(reflection, -surface->ray->direction), 0, 1);

    phong_color +=
        PhongColor(lights[l].color / (float)(config.extra_lights_per_light + 1),
                   color_diffuse, color_specular, l_dot_n, r_dot_v, shininess);

    // launch extra shadow rays
    for (uint e = 0; e < config.extra_lights_per_light; ++e) {
      Ray shadow_ray;
      shadow_ray.position = intersection_position + kShadowRayBias * normal;

      uint k = l * config.extra_lights_per_light + e;

      glm::vec3 shadow_to_light =
          extra_lights->positions[k] - shadow_ray.position;
      float shadow_to_light_dist = glm::length(shadow_to_light);

      shadow_ray.direction = shadow_to_light / shadow_to_light_dist;

      int hit = Intersect(&shadow_ray, spheres, sphere_count, triangles,
                          triangle_count, &occluder);

      int occluded = hit && occluder.t + kTolerance < shadow_to_light_dist;
      if (occluded) {
        continue;
      }

      l_dot_n = glm::clamp<float>(glm::dot(shadow_ray.direction, normal), 0, 1);
      reflection = glm::normalize(2 * l_dot_n * normal - shadow_ray.direction);
      r_dot_v = glm::clamp<float>(
          glm::dot(reflection, -surface->ray->direction), 0, 1);

      phong_color += PhongColor(
          lights[l].color / (float)(config.extra_lights_per_light + 1),
          color_diffuse, color_specular, l_dot_n, r_dot_v, shininess);
    }
  }

  // Launch reflected rays.
  if (bounces > 0) {
    Ray ray;
    ray.position = intersection_position;
    float vn =
        glm::clamp<float>(glm::dot(-surface->ray->direction, normal), 0, 1);
    ray.direction = glm::normalize(2 * vn * normal + surface->ray->direction);
    ray.position += kReflectionRayBias * normal;

    glm::vec3 color = TraceRay(&ray, scene, bounces - 1, extra_lights);
    for (int i = 0; i < kRgbChannel__Count; ++i) {
      phong_color[i] *= (1 - color_specular[i]);
      phong_color[i] += color_specular[i] * color[i];
    }
  }

  return phong_color;
}

static Status ParseConfig(uint argc, char *argv[], Config *c) {
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

  Scene scene;
  st = LoadScene(config.scene_filepath, &ambient_light, &scene);
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
  glutIdleFunc(Idle);
  Init();

  Lights extra_lights;

  if (config.extra_lights_per_light > 0) {
    std::default_random_engine re(std::random_device{}());
    std::uniform_real_distribution<float> distrib(-1, 1);

    uint count = scene.light_count * config.extra_lights_per_light;
    extra_lights.positions.resize(count);
    extra_lights.colors.resize(count);

    for (uint i = 0; i < scene.light_count; ++i) {
      for (uint j = 0; j < config.extra_lights_per_light; ++j) {
        float x = distrib(re);
        float y = distrib(re);
        glm::vec3 offset(x, y, 0);

        uint k = i * config.extra_lights_per_light + j;

        extra_lights.positions[k] = scene.lights[i].position + offset;
        extra_lights.colors[k] = scene.lights[i].color;
      }
    }
  }

  std::vector<Ray> rays;
  if (config.rays_per_pixel == 1) {
    MakeDefaultRays(&kCameraPosition, IMG_W, IMG_H, FOV, FOCAL_LEN, &rays);
  } else if (config.rays_per_pixel > 1) {
    MakeJitteredRays(&kCameraPosition, IMG_W, IMG_H, FOV, FOCAL_LEN,
                     config.rays_per_pixel, &rays);
  }

  for (int y = 0; y < IMG_H; ++y) {
    for (int x = 0; x < IMG_W; ++x) {
      glm::vec3 color;
      for (int i = 0; i < config.rays_per_pixel; ++i) {
        uint j = (y * IMG_W + x) * config.rays_per_pixel + i;
        color += TraceRay(&rays[j], &scene, config.reflection_bounces,
                          &extra_lights);
      }
      for (int i = 0; i < kRgbChannel__Count; ++i) {
        buffer[y][x][i] =
            glm::clamp<float>(color[i] / config.rays_per_pixel, 0, 1) * 255;
      }
    }
  }

  glutMainLoop();
}