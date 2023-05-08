#ifndef TRAYCER_MAIN_HPP
#define TRAYCER_MAIN_HPP

#include <glm/glm/vec3.hpp>
#include <vector>

#define FILEPATH_BUFFER_SIZE 4096

enum Status {
  kStatus_Ok,
  kStatus_UnspecifiedError,
  kStatus_IoError,
  kStatus_GlError
};

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

struct Vertex {
  glm::vec3 position;
  glm::vec3 color_diffuse;
  glm::vec3 color_specular;
  glm::vec3 normal;
  float shininess;
};

struct Triangle {
  Vertex v[3];
};

struct Sphere {
  glm::vec3 position;
  glm::vec3 color_diffuse;
  glm::vec3 color_specular;
  float shininess;
  float radius;
};

struct Light {
  glm::vec3 position;
  glm::vec3 color;
};

struct Lights {
  std::vector<glm::vec3> positions;
  std::vector<glm::vec3> colors;
};

struct Ray {
  glm::vec3 position;
  glm::vec3 direction;
};

enum ObjectType { SPHERE, TRIANGLE };

struct TriangleData {
  int index;
  float alpha;
  float beta;
  float gamma;
};

struct SphereData {
  int index;
};

union ObjectData {
  SphereData sphere;
  TriangleData triangle;
};

struct Intersection {
  ObjectType type;
  ObjectData data;
  Ray *ray;
  float t;
  bool hit;
};

enum SamplingMode { DEFAULT, SUPER_JITTER };

enum RenderMode { kRenderMode_Display, kRenderMode_Jpeg, kRenderMode__Count };

#endif  // TRAYCER_MAIN_HPP