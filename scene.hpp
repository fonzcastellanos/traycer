#ifndef TRAYCER_SCENE_HPP
#define TRAYCER_SCENE_HPP

#include <glm/glm/vec3.hpp>

#include "status.hpp"
#include "types.hpp"

#define MAX_LIGHT_COUNT 128
#define MAX_SPHERE_COUNT 128
#define MAX_TRIANGLE_COUNT 32768

struct Sphere {
  glm::vec3 position;
  glm::vec3 color_diffuse;
  glm::vec3 color_specular;
  float shininess;
  float radius;
};

struct Vertex {
  glm::vec3 position;
  glm::vec3 color_diffuse;
  glm::vec3 color_specular;
  glm::vec3 normal;
  float shininess;
};

struct Triangle {
  Vertex vertices[3];
  glm::vec3 normal;
  float area;
};

struct TriangleVertices {
  glm::vec3 position[MAX_TRIANGLE_COUNT * 3];
  glm::vec3 diffuse_color[MAX_TRIANGLE_COUNT * 3];
  glm::vec3 specular_color[MAX_TRIANGLE_COUNT * 3];
  glm::vec3 normal[MAX_TRIANGLE_COUNT * 3];
  float shininess[MAX_TRIANGLE_COUNT * 3];
};

struct Triangles {
  TriangleVertices vertices;
  glm::vec3 normal[MAX_TRIANGLE_COUNT];
  float area[MAX_TRIANGLE_COUNT];
  uint count;
};

struct Light {
  glm::vec3 position;
  glm::vec3 color;
};

struct Scene {
  glm::vec3 ambient_light;

  Light lights[MAX_LIGHT_COUNT];
  uint light_count;

  Sphere spheres[MAX_SPHERE_COUNT];
  uint sphere_count;

  Triangles triangles;
};

Status LoadScene(const char* filepath, Scene* scene);

#endif  // TRAYCER_SCENE_HPP