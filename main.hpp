#ifndef TRAYCER_MAIN_HPP
#define TRAYCER_MAIN_HPP

#include <glm/glm/vec3.hpp>

enum Status {
  kStatus_Ok,
  kStatus_UnspecifiedError,
  kStatus_IoError,
  kStatus_GlError
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