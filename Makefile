RTRACE_CXX_SRC=rtrace.cpp
RTRACE_HEADER=
RTRACE_OBJ=$(notdir $(patsubst %.cpp,%.o,$(RTRACE_CXX_SRC)))

IMAGE_LIB_SRC=$(wildcard external/imageIO/*.cpp)
IMAGE_LIB_HEADER=$(wildcard external/imageIO/*.h)
IMAGE_LIB_OBJ=$(notdir $(patsubst %.cpp,%.o,$(IMAGE_LIB_SRC)))

HEADER=$(RTRACE_HEADER) $(IMAGE_LIB_HEADER)
CXX_OBJ=$(RTRACE_OBJ) $(IMAGE_LIB_OBJ)

CXX=g++
TARGET=rtrace
CXXFLAGS=-DGLM_FORCE_RADIANS -Wno-unused-result
OPT=-O3

UNAME_S=$(shell uname -s)

ifeq ($(UNAME_S),Linux)
  PLATFORM=Linux
  INCLUDE=-Iexternal/glm/ -Iexternal/imageIO
  LIB=-lGLEW -lGL -lglut -ljpeg
  LDFLAGS=
else
  PLATFORM=Mac OS
  INCLUDE=-Iexternal/glm/ -Iexternal/imageIO -Iexternal/jpeg-9a-mac/include
  LIB=-framework OpenGL -framework GLUT external/jpeg-9a-mac/lib/libjpeg.a
  CXXFLAGS+= -Wno-deprecated-declarations
  LDFLAGS=-Wl,-w
endif

all: $(TARGET)

$(TARGET): $(CXX_OBJ)
	$(CXX) $(LDFLAGS) $^ $(OPT) $(LIB) -o $@

$(RTRACE_OBJ):%.o: %.cpp $(HEADER)
	$(CXX) -c $(CXXFLAGS) $(OPT) $(INCLUDE) $< -o $@

$(IMAGE_LIB_OBJ):%.o: external/imageIO/%.cpp $(IMAGE_LIB_HEADER)
	$(CXX) -c $(CXXFLAGS) $(OPT) $(INCLUDE) $< -o $@

clean:
	rm -rf *.o $(TARGET)
