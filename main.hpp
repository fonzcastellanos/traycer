#ifndef TRAYCER_MAIN_HPP
#define TRAYCER_MAIN_HPP

#include <glm/glm/vec3.hpp>
#include <vector>

#define FILEPATH_BUFFER_SIZE 4096

struct Config {
  char *scene_filepath;
  int rays_per_pixel;
  int reflection_bounces;
  uint extra_lights_per_light;
  char render_filepath[FILEPATH_BUFFER_SIZE];
};

enum RgbChannel {
  kRgbChannel_Red,
  kRgbChannel_Green,
  kRgbChannel_Blue,
  kRgbChannel__Count
};

struct Lights {
  std::vector<glm::vec3> positions;
  std::vector<glm::vec3> colors;
};

struct Ray {
  glm::vec3 position;
  glm::vec3 direction;
};

enum GeometryType { kGeometryType_Sphere, kGeometryType_Triangle };

struct TriangleIntersection {
  uint index;
  float alpha;
  float beta;
  float gamma;
};

struct SphereIntersection {
  uint index;
};

struct Intersection {
  GeometryType type;
  union {
    SphereIntersection sphere;
    TriangleIntersection triangle;
  };
  Ray *ray;
  float t;
  bool hit;
};

enum SamplingMode { kSamplingMode_Default, kSamplingMode_Jitter };

enum RenderTarget {
  kRenderTarget_Window,
  kRenderTarget_Jpeg,
  kRenderTarget__Count
};

#endif  // TRAYCER_MAIN_HPP