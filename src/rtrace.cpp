/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: <Alfonso Castellanos>
 * *************************
*/

#ifdef WIN32
    #include <windows.h>
#endif

#if defined(WIN32) || defined(linux)
    #include <GL/gl.h>
    #include <GL/glut.h>
#elif defined(__APPLE__)
    #include <OpenGL/gl.h>
    #include <GLUT/glut.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef WIN32
    #define strcasecmp _stricmp
#endif

#include <imageIO.h>
#include <glm/glm.hpp>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <cfloat>

#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100

char * filename = NULL;

// different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode = MODE_DISPLAY;

// you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

const int numPixels = WIDTH*HEIGHT;
const float aspectRatio = (float)WIDTH/HEIGHT;

const float eps = 0.00000001;
const glm::vec3 backgroundColor(1.0, 1.0, 1.0);

int raysPerPixel = 0;
int reflectionBounces = 0;

// camera field of view
#define fov 60.0

// camera position
glm::vec3 camPos(0.0, 0.0, 0.0);

unsigned char buffer[HEIGHT][WIDTH][3];

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

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
glm::vec3 ambient_light;

int numExtraLights;
Light * extraLights[MAX_LIGHTS];

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;

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

union ObjectData { SphereData sphere; TriangleData triangle; };

struct Intersection {
    ObjectType type;
    ObjectData data;
    Ray* ray;
    float t;
    bool hit;
};

enum SamplingMode { DEFAULT, SUPER_JITTER };

void draw_scene();
void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void save_jpg();
void parse_check(const char *expected, char *found);
void parse_floats(FILE* file, const char *check, glm::vec3 & p);
void parse_rad(FILE *file, float *r);
void parse_shi(FILE *file, float *shi);
int loadScene(char *argv);
void display();
void init();
void idle();
Ray ** makeRays(SamplingMode mode);
void getReflectedRay(Ray * ray, Intersection * in, Ray * reflectedRay);
void sphereIntersect(Ray * ray, int s, Intersection * out);
void triangleIntersect(Ray * ray, int tri, Intersection * out);
void intersect(Ray * ray, Intersection * prev, Intersection * out);
float getPhongColor(float lightColor, float diffuse, float specular, float ln, float rv, float shininess);
glm::vec3 traceRay(Ray * ray, Intersection * prev, int bounces);
glm::vec3 shade(Intersection * surface, int bounces);

void draw_scene() {
    for(unsigned int x=0; x<WIDTH; x++) {
        glPointSize(2.0);  
        glBegin(GL_POINTS);
        for(unsigned int y=0; y<HEIGHT; y++) {
            plot_pixel(x, y, buffer[y][x][0], buffer[y][x][1], buffer[y][x][2]);
        }
        glEnd();
        glFlush();
    }
    printf("Done!\n"); fflush(stdout);
}

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b) {
    glColor3f(((float)r) / 255.0f, ((float)g) / 255.0f, ((float)b) / 255.0f);
    glVertex2i(x,y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b) {
    buffer[y][x][0] = r;
    buffer[y][x][1] = g;
    buffer[y][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b) {
    plot_pixel_display(x,y,r,g,b);
    if(mode == MODE_JPEG)
        plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg() {
    printf("Saving JPEG file: %s\n", filename);

    ImageIO img(WIDTH, HEIGHT, 3, &buffer[0][0][0]);
    if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
        printf("Error in Saving\n");
    else 
        printf("File saved Successfully\n");
}

void parse_check(const char *expected, char *found) {
    if(strcasecmp(expected,found)) {
        printf("Expected '%s ' found '%s '\n", expected, found);
        printf("Parse error, abnormal abortion\n");
        exit(0);
    }
}

void parse_floats(FILE* file, const char *check, glm::vec3 & p) {
    char str[100];
    fscanf(file,"%s",str);
    parse_check(check,str);
    fscanf(file,"%f %f %f",&p[0],&p[1],&p[2]);
    printf("%s %f %f %f\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE *file, float *r) {
    char str[100];
    fscanf(file,"%s",str);
    parse_check("rad:",str);
    fscanf(file,"%f",r);
    printf("rad: %f\n",*r);
}

void parse_shi(FILE *file, float *shi) {
    char s[100];
    fscanf(file,"%s",s);
    parse_check("shi:",s);
    fscanf(file,"%f",shi);
    printf("shi: %f\n",*shi);
}

int loadScene(char *argv) {
    FILE * file = fopen(argv,"r");
    int number_of_objects;
    char type[50];
    Triangle t;
    Sphere s;
    Light l;
    fscanf(file,"%i", &number_of_objects);

    printf("number of objects: %i\n",number_of_objects);

    parse_floats(file,"amb:",ambient_light);

    for(int i=0; i<number_of_objects; i++) {
        fscanf(file,"%s\n",type);
        printf("%s\n",type);
        if(strcasecmp(type,"triangle")==0) {
            printf("found triangle\n");
            for(int j=0;j < 3;j++) {
                parse_floats(file,"pos:",t.v[j].position);
                parse_floats(file,"nor:",t.v[j].normal);
                parse_floats(file,"dif:",t.v[j].color_diffuse);
                parse_floats(file,"spe:",t.v[j].color_specular);
                parse_shi(file,&t.v[j].shininess);
            }

            if(num_triangles == MAX_TRIANGLES) {
                printf("too many triangles, you should increase MAX_TRIANGLES!\n");
                exit(0);
            }
            triangles[num_triangles++] = t;
        }
        else if(strcasecmp(type,"sphere")==0) {
            printf("found sphere\n");

            parse_floats(file,"pos:",s.position);
            parse_rad(file,&s.radius);
            parse_floats(file,"dif:",s.color_diffuse);
            parse_floats(file,"spe:",s.color_specular);
            parse_shi(file,&s.shininess);

            if(num_spheres == MAX_SPHERES) {
                printf("too many spheres, you should increase MAX_SPHERES!\n");
                exit(0);
            }
            spheres[num_spheres++] = s;
        }
        else if(strcasecmp(type,"light")==0) {
            printf("found light\n");
            parse_floats(file,"pos:",l.position);
            parse_floats(file,"col:",l.color);

            if(num_lights == MAX_LIGHTS) {
                printf("too many lights, you should increase MAX_LIGHTS!\n");
                exit(0);
            }
            lights[num_lights++] = l;
        }
        else {
            printf("unknown type in scene description:\n%s\n",type);
            exit(0);
        }
    }
    return 0;
}

void display() {
}

void init() {
    glMatrixMode(GL_PROJECTION);
    glOrtho(0,WIDTH,0,HEIGHT,1,-1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glClearColor(0,0,0,0);
    glClear(GL_COLOR_BUFFER_BIT);
}

void idle() {
    //hack to make it only draw once
    static int once=0;
    if(!once) {
        draw_scene();
        if(mode == MODE_JPEG)
            save_jpg();
    }
    once=1;
}

Ray ** makeRays(SamplingMode mode) {
    // Generate image plane corners
    const float z = -1.0;
    float tanHalfFOV = glm::tan( (fov/2.0) * (M_PI/180.0) );
    glm::vec3 topRight(aspectRatio*tanHalfFOV, tanHalfFOV, z);
    glm::vec3 topLeft(-aspectRatio*tanHalfFOV, tanHalfFOV, z);
    glm::vec3 bottomRight(aspectRatio*tanHalfFOV, -tanHalfFOV, z);
    glm::vec3 bottomLeft(-aspectRatio*tanHalfFOV, -tanHalfFOV, z);

    // Allocate ray array
    Ray ** rays = new Ray * [numPixels];
    for (int i = 0; i < numPixels; ++i) {
        rays[i] = new Ray[raysPerPixel];
    }

    // Calculate pixel dimensions
    float pWidth = (topRight.x - topLeft.x)/WIDTH;
    float pHeight = (topRight.y - bottomRight.y)/HEIGHT;

    // Sample
    int rayInd = 0;
    if (mode == DEFAULT) {
        for (float y = bottomRight.y; y < topRight.y-pHeight/2.0; y += pHeight) {
            for (float x = topLeft.x; x < topRight.x-pWidth/2.0; x += pWidth) {
                rays[rayInd][0].position = glm::vec3(x+pWidth/2.0, y+pHeight/2.0, z);
                rays[rayInd][0].direction = glm::normalize(rays[rayInd][0].position - camPos);
                ++rayInd;
            }
        }
    }
    else if (mode == SUPER_JITTER) {
        std::srand(std::time(NULL));
        for (float y = bottomRight.y; y < topRight.y-pHeight/2.0; y += pHeight) {
            for (float x = topLeft.x; x < topRight.x-pWidth/2.0; x += pWidth) {
                for (int r = 0; r < raysPerPixel; ++r) {
                    float xOffset = pWidth*((float)(std::rand())/(float)(RAND_MAX)); 
                    float yOffset = pHeight*((float)(std::rand())/(float)(RAND_MAX));
                    rays[rayInd][r].position = glm::vec3(x+xOffset, y+yOffset, z);
                    rays[rayInd][r].direction = glm::normalize(rays[rayInd][r].position - camPos);
                }
                ++rayInd;
            }
        }
    }
    return rays;
}

void getReflectedRay(Ray * ray, Intersection * in, Ray * reflectedRay) {
    glm::vec3 p = ray->position + in->t*ray->direction;
    glm::vec3 n;
    if (in->type == SPHERE) {
        int s = in->data.sphere.index;
        n = (p-spheres[s].position)/spheres[s].radius; 
    }
    else if (in->type == TRIANGLE) {
        int tri = in->data.triangle.index;
        float alpha = in->data.triangle.alpha;
        float beta = in->data.triangle.beta;
        float gamma = in->data.triangle.gamma;

        n = triangles[tri].v[0].normal*alpha+
            triangles[tri].v[1].normal*beta+
            triangles[tri].v[2].normal*gamma;
        n = glm::normalize(n);

    }
    float vn = glm::clamp<float>(glm::dot(-ray->direction, n), 0.0, 1.0);
    reflectedRay->direction = glm::normalize(2.0f*vn*n+ray->direction);
    reflectedRay->position = p;
}

void sphereIntersect(Ray * ray, int s, Intersection * out) {
    if (ray == NULL || out == NULL) return;

    out->type = SPHERE;
    out->data.sphere.index = s;
    out->ray = ray;
    out->hit = false;

    // Calculate quadratic coefficients
    glm::vec3 centerToRay = ray->position - spheres[s].position;
    float b = 2.0*glm::dot(ray->direction, centerToRay);
    float c = glm::pow(centerToRay.x, 2) + glm::pow(centerToRay.y, 2) + glm::pow(centerToRay.z, 2) - glm::pow(spheres[s].radius, 2);

    // Calculate discriminant
    float discriminant = glm::pow(b, 2)-4.0*c;

    // Calculate roots and determine t
    float t;
    if (discriminant < -eps) { // misses sphere
        return;
    }
    else if (glm::abs(discriminant) < eps) {  // hits sphere at one point
        t = -b/2.0;
    }
    else { // hits sphere at two points
        float q = -0.5*(b+glm::sign(b)*glm::sqrt(discriminant));
        float root1 = c/q;
        float root2 = q;
        t = glm::min(root1, root2);
    }

    if (t > eps) { // hit because ahead of origin
        out->t = t;
        out->hit = true;
    }
}

void triangleIntersect(Ray * ray, int tri, Intersection * out) {
    if (ray == NULL || out == NULL) return;

    out->type = TRIANGLE;
    out->data.triangle.index = tri;
    out->ray = ray;
    out->hit = false;

    // Alias for vertices
    Vertex (&v)[3] = triangles[tri].v;

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
    float t = -(glm::dot(n, ray->position)+d)/nDotRay;
    if (t < eps) return;

    glm::vec3 p = ray->position + t*ray->direction;

    glm::vec3 np; // vector perpendicular to subtriangle

    // Inside-outside test
    glm::vec3 edge0p = p - v[0].position;
    np = glm::cross(edge01, edge0p);
    if (glm::dot(n, np) < -eps) return;

    float gamma = glm::length(np)/twiceTriArea;

    glm::vec3 edge12 = v[2].position - v[1].position;
    glm::vec3 edge1p = p - v[1].position;
    np = glm::cross(edge12, edge1p);
    if (glm::dot(n, np) < -eps) return;

    float alpha = glm::length(np)/twiceTriArea; 

    glm::vec3 edge20 = v[0].position - v[2].position;
    glm::vec3 edge2p = p - v[2].position;
    np = glm::cross(edge20, edge2p);
    if (glm::dot(n, np) < -eps) return;

    float beta = glm::length(np)/twiceTriArea; 

    out->data.triangle.alpha = alpha;
    out->data.triangle.beta = beta;
    out->data.triangle.gamma = gamma;
    out->t = t;
    out->hit = true;
}

void intersect(Ray * ray, Intersection * prev, Intersection * out) {
    if (ray == NULL || out == NULL) return;
    out->ray = ray;
    out->t = (float)FLT_MAX;
    out->hit = false;

    Intersection current;
    // sphere intersections
    for (int s = 0; s < num_spheres; ++s) {
        if (prev != NULL && prev->type == SPHERE && prev->data.sphere.index == s) continue;

        sphereIntersect(ray, s, &current);

        if (!current.hit || current.t > out->t+eps) continue;

        out->type = current.type;
        out->data.sphere.index = current.data.sphere.index;
        out->ray = current.ray;
        out->t = current.t;
        out->hit = current.hit;
    }
    // triangle intersections
    for (int tri = 0; tri < num_triangles; ++tri) {
        if (prev != NULL && prev->type == TRIANGLE && prev->data.triangle.index == tri) continue;

        triangleIntersect(ray, tri, &current);

        if (!current.hit || current.t > out->t+eps) continue;

        out->type = current.type;
        out->data.triangle.index = current.data.triangle.index;
        out->data.triangle.alpha= current.data.triangle.alpha;
        out->data.triangle.beta = current.data.triangle.beta;
        out->data.triangle.gamma = current.data.triangle.gamma;
        out->ray = current.ray;
        out->t = current.t;
        out->hit = current.hit;
    }
}

float getPhongColor(float lightColor, float diffuse, float specular, float ln, float rv, float shininess) {
    return lightColor*(diffuse*ln + specular*glm::pow(rv, shininess));
}


glm::vec3 traceRay(Ray * ray, Intersection * prev, int bounces) {
    Intersection closest;
    intersect(ray, prev, &closest); 
    if (!closest.hit) {
        return backgroundColor; 
    }
    else {
        return shade(&closest, bounces);
    }
}

glm::vec3 shade(Intersection * surface, int bounces) {
    int index;
    glm::vec3 p; // point of intersection
    glm::vec3 n; // normal;
    glm::vec3 color_diffuse;
    glm::vec3 color_specular;
    float shininess;
    
    // calculate surface-specific parameters
    if (surface->type == SPHERE) {
        index = surface->data.sphere.index;
        p = surface->ray->position + surface->t*(surface->ray->direction); // point of intersection
        n = (p-spheres[index].position)/spheres[index].radius;  // normal
        color_diffuse = spheres[index].color_diffuse;
        color_specular = spheres[index].color_specular;
        shininess = spheres[index].shininess;
    }
    else if (surface->type == TRIANGLE) {
        index = surface->data.triangle.index;
        p = surface->ray->position + surface->t*(surface->ray->direction);

        float alpha = surface->data.triangle.alpha;
        float beta = surface->data.triangle.beta;
        float gamma = surface->data.triangle.gamma;

        n = triangles[index].v[0].normal*alpha+
            triangles[index].v[1].normal*beta+
            triangles[index].v[2].normal*gamma;
        n = glm::normalize(n);

        color_diffuse = triangles[index].v[0].color_diffuse*alpha+
                        triangles[index].v[1].color_diffuse*beta+
                        triangles[index].v[2].color_diffuse*gamma;

        color_specular = triangles[index].v[0].color_specular*alpha+
                         triangles[index].v[1].color_specular*beta+
                         triangles[index].v[2].color_specular*gamma;

        shininess = triangles[index].v[0].shininess*alpha+
                    triangles[index].v[1].shininess*beta+
                    triangles[index].v[2].shininess*gamma;
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
        shadow.direction = glm::normalize(lights[l].position-p);

        intersect(&shadow, surface, &occluder);

        toLight = glm::length(lights[l].position-shadow.position);
        if (occluder.hit && occluder.t+eps < toLight) continue;

        ln = glm::clamp<float>(glm::dot(shadow.direction, n), 0.0, 1.0);
        reflection = glm::normalize(2.0f*ln*n-shadow.direction);
        rv = glm::clamp<float>(glm::dot(reflection, -surface->ray->direction), 0.0, 1.0);

        for (int c = 0; c < 3; ++c) {
            phongColor[c] += getPhongColor(lights[l].color[c]/(numExtraLights+1), 
                                           color_diffuse[c], 
                                           color_specular[c], 
                                           ln, rv, shininess);
        }

        // launch extra shadow rays
        for (int e = 0; e < numExtraLights; ++e) {
            shadow.position = p;
            shadow.direction = glm::normalize(extraLights[l][e].position-p);

            intersect(&shadow, surface, &occluder);

            toLight = glm::length(extraLights[l][e].position-shadow.position);
            if (occluder.hit && occluder.t+eps < toLight) continue;

            ln = glm::clamp<float>(glm::dot(shadow.direction, n), 0.0, 1.0);
            reflection = glm::normalize(2.0f*ln*n-shadow.direction);
            rv = glm::clamp<float>(glm::dot(reflection, -surface->ray->direction), 0.0, 1.0);

            for (int c = 0; c < 3; ++c) {
                phongColor[c] += getPhongColor(lights[l].color[c]/(numExtraLights+1), 
                                               color_diffuse[c], 
                                               color_specular[c], 
                                               ln, rv, shininess);
            }
        }
    }

    // launch reflected rays
    if (bounces > 0) {
        Ray reflectedRay;
        getReflectedRay(surface->ray, surface, &reflectedRay);
        glm::vec3 reflectedColor = traceRay(&reflectedRay, surface, bounces-1);
        for (int c = 0; c < 3; ++c) {
            phongColor[c] *= (1-color_specular[c]);
            phongColor[c] += color_specular[c]*reflectedColor[c];
        }
    }

    return phongColor;
}

int main(int argc, char ** argv) {
    if ((argc < 5) || (argc > 6)) {  
    printf ("Usage: %s <input scenefile> <rays per pixel> <reflection bounces> <extra lights per light> [output jpegname] \n", argv[0]);
    exit(0);
    }
    if(argc == 6) {
        mode = MODE_JPEG;
        filename = argv[5];
    }
    else if(argc == 5) {
        mode = MODE_DISPLAY;
    }
    raysPerPixel = std::atoi(argv[2]);
    reflectionBounces = std::atoi(argv[3]);
    numExtraLights = std::atoi(argv[4]);

    glutInit(&argc,argv);
    loadScene(argv[1]);

    glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
    glutInitWindowPosition(0,0);
    glutInitWindowSize(WIDTH,HEIGHT);
    int window = glutCreateWindow("Ray Tracer");
    glutDisplayFunc(display);
    glutIdleFunc(idle);
    init();

    if (numExtraLights > 0) {
        std::srand(std::time(NULL));
        for (int l = 0; l < num_lights; ++l) {
            extraLights[l] = new Light[numExtraLights];
            for (int e = 0; e < numExtraLights; ++e) {
                float xOffset = (-1.0 + 2.0*((float)(std::rand())/(float)(RAND_MAX))); 
                float yOffset = (-1.0 + 2.0*((float)(std::rand())/(float)(RAND_MAX)));
                //float zOffset = (-1.0 + 2.0*((float)(std::rand())/(float)(RAND_MAX)));
                glm::vec3 offset(xOffset, yOffset, 0.0);

                extraLights[l][e].position = lights[l].position + offset;
                extraLights[l][e].color = lights[l].color;
            }
        }
    }

    Ray ** rays = NULL;
    if (raysPerPixel == 1) {
        rays = makeRays(DEFAULT);
    }
    else if (raysPerPixel > 1) {
        rays = makeRays(SUPER_JITTER);
    }

    for (int p = 0; p < numPixels; ++p) {
        int x = p%WIDTH;
        int y = p/WIDTH;
        glm::vec3 totalColor(0.0, 0.0, 0.0);
        for (int r = 0; r < raysPerPixel; ++r) {
            totalColor += traceRay(&rays[p][r], NULL, reflectionBounces);
        }
        for (int c = 0; c < 3; ++c) {
            buffer[y][x][c] = glm::clamp<float>(totalColor[c]/raysPerPixel, 0.0, 1.0)*255.0;
        }
    }
    delete [] rays;

    glutMainLoop();
}

