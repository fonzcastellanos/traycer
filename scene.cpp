#include "scene.hpp"

#include <cassert>
#include <cstdio>
#include <cstring>

Status ParseField(std::FILE *f, const char *name, glm::vec3 *vals) {
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

Status LoadScene(const char *filepath, glm::vec3 *ambient_light,
                 std::vector<Triangle> *triangles, std::vector<Sphere> *spheres,
                 Light *lights, uint *light_count) {
  assert(filepath);
  assert(triangles);
  assert(spheres);
  assert(lights);
  assert(light_count);

  std::FILE *file = std::fopen(filepath, "r");
  if (!file) {
    std::fprintf(stderr, "Failed to open file %s.\n", filepath);
    return kStatus_IoError;
  }

  uint obj_count;
  char type[50];
  Light l;

  int rc = std::fscanf(file, "%u", &obj_count);
  if (rc != 1) {
    std::fprintf(stderr, "Failed to parse object count from file %s.\n",
                 filepath);
    return kStatus_IoError;
  }

  std::printf("Object count: %u\n", obj_count);

  Status st = ParseField(file, "amb:", ambient_light);
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
      Triangle t;
      st = ParseTriangle(file, &t);
      if (st != kStatus_Ok) {
        std::fprintf(stderr, "Failed to parse triangle.\n");
        return st;
      }

      triangles->push_back(t);
    } else if (std::strcmp(type, "sphere") == 0) {
      Sphere sph;
      st = ParseSphere(file, &sph);
      if (st != kStatus_Ok) {
        std::fprintf(stderr, "Failed to parse sphere.\n");
        return st;
      }

      spheres->push_back(sph);
    } else if (std::strcmp(type, "light") == 0) {
      st = ParseLight(file, &l);
      if (st != kStatus_Ok) {
        std::fprintf(stderr, "Failed to parse light.\n");
        return st;
      }

      if (*light_count == MAX_LIGHTS) {
        std::fprintf(stderr,
                     "Too many lights. Increase MAX_LIGHTS if more lights are "
                     "desired.\n");
        // TODO: Choose a better-suited status.
        return kStatus_UnspecifiedError;
      }

      lights[*light_count] = l;
      ++(*light_count);
    } else {
      std::fprintf(stderr, "Invalid object type \"%s\" in scene file %s\n",
                   type, filepath);

      // TODO: Choose a better-suited status.
      return kStatus_UnspecifiedError;
    }
  }

  return kStatus_Ok;
}
