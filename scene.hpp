#ifndef TRAYCER_SCENE_HPP
#define TRAYCER_SCENE_HPP

#include <glm/glm/vec3.hpp>
#include <vector>

#include "status.hpp"
#include "types.hpp"

#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100

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
  Vertex v[3];
};

struct Light {
  glm::vec3 position;
  glm::vec3 color;
};

Status LoadScene(const char* filepath, glm::vec3* ambient_light,
                 std::vector<Triangle>* triangles, std::vector<Sphere>* spheres, Light* lights, uint* light_count);

#endif  // TRAYCER_SCENE_HPP