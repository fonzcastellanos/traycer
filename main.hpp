#ifndef TRAYCER_MAIN_HPP
#define TRAYCER_MAIN_HPP

#include <glm/glm/vec3.hpp>
#include <vector>

#define FILEPATH_BUFFER_SIZE 4096

struct Config {
  char *scene_filepath;
  int jitter;
  int bounces;
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

enum SurfaceType {
  kSurfaceType_Sphere,
  kSurfaceType_Triangle,
  kSurfaceType__Count
};

struct TriangleIntersection {
  float alpha;
  float beta;
  float gamma;
};

struct Intersection {
  SurfaceType type;
  uint index;
  float t;
  TriangleIntersection triangle;
};

enum SamplingMode { kSamplingMode_Default, kSamplingMode_Jitter };

enum RenderTarget {
  kRenderTarget_Window,
  kRenderTarget_Jpeg,
  kRenderTarget__Count
};

#endif  // TRAYCER_MAIN_HPP