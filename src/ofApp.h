//
//  RayCaster - Set of simple classes to create a camera/view setup for our Ray Tracer HW Project
//
//  I've included these classes as a mini-framework for our introductory ray tracer. The
//  classes are written in simplified C++ for those who are transitioning from Java.
//  You are free to modify/change.   
//
//  These classes provide a simple render camera which can can return a ray starting from
//  it's position to a (u, v) coordinate on the view plane.
//
//  The view plane is where we can locate our photorealistic image we are rendering.
//  The field-of-view of the camera by moving it closer/further 
//  from the view plane.  The viewplane can be also resized.  When ray tracing an image, the aspect
//  ratio of the view plane should the be same as your image. So for example, the current view plane
//  default size is ( 6.0 width by 4.0 height ).   A 1200x800 pixel image would have the same
//  aspect ratio.
//
//  This is not a complete ray tracer - just a set of skelton classes to start.  The current
//  base scene object only stores a value for the diffuse/specular color of the object (defaut is gray).
//  at some point, we will want to replace this with a Material class that contains these (and other 
//  parameters)
//  
//  (c) Kevin M. Smith  - 24 September 2018
//  Calvin Quach - 21 May 2023
//  - created additional classes for Toruses, boxes, and 3D fractals
//  - implemented sdf's and other needed parameters for each new object
#pragma once

#include "ofMain.h"
#include <glm/gtx/intersect.hpp>
#include <glm/gtc/noise.hpp>


#define MAX_RAY_STEPS 1000
#define DIST_THRESHOLD .01
#define MAX_DISTANCE 100
#define RAY_EPS .15

//  General Purpose Ray class 
//
class Ray {
public:
	Ray(glm::vec3 p, glm::vec3 d) { this->p = p; this->d = d; }
	void draw(float t) { ofDrawLine(p, p + t * d); }

	glm::vec3 evalPoint(float t) {
		return (p + t * d);
	}

	glm::vec3 p, d;
};

//  Base class for any renderable object in the scene
//
class SceneObject {
public: 
	virtual void draw() = 0;    // pure virtual funcs - must be overloaded
	virtual bool intersect(const Ray &ray, glm::vec3 &point, glm::vec3 &normal) { return false; }
	virtual float sdf(const glm::vec3 &point) { return 0.0; }
	virtual glm::vec3 getNormal(const glm::vec3 &p) { return glm::vec3(1, 0, 0); }
	
	// any data common to all scene objects goes here
	glm::vec3 position = glm::vec3(0, 0, 0);

	// material properties (we will ultimately replace this with a Material class - TBD)
	//
	ofColor diffuseColor = ofColor::grey;    // default colors - can be changed.
	ofColor specularColor = ofColor::lightGray;
	bool isSelectable = true;
};

//  General purpose sphere  (assume parametric)
//
class Sphere: public SceneObject {
public:
	Sphere(glm::vec3 p, float r, ofColor diffuse = ofColor::lightGray) { position = p; radius = r; diffuseColor = diffuse; }
	Sphere() {}
	bool intersect(const Ray &ray, glm::vec3 &point, glm::vec3 &normal) {
		return (glm::intersectRaySphere(ray.p, ray.d, position, radius, point, normal));
	}

	float sdf(const glm::vec3 &p) { 
		return glm::length(position - p) - radius; 
	}

	glm::vec3 getNormal(const glm::vec3 &p) { return glm::normalize(p - position); }

	void draw()  { 
	//	ofDrawSphere(position, radius); 
		spherePrim.setRadius(radius);
		spherePrim.setPosition(position);
		spherePrim.draw();
	}

	float radius = 1.0;
	ofSpherePrimitive spherePrim;
};

/*
* Primitive Torus Class
*/
class Torus : public SceneObject {
public:
	Torus(glm::vec3 p, ofColor diffuse = ofColor::lightGray) { 
		position = p; 
		diffuseColor = diffuse; 

		trans = glm::translate(glm::mat4(1.0), position);
		rot1 = glm::rotate(trans, glm::radians(90.0f), glm::vec3(1, 0, 0));
		rot2 = glm::rotate(rot1, glm::radians(0.0f), glm::vec3(0, 1, 0));
		sca = glm::scale(rot2, glm::vec3(1, 1, 1));
	}
	Torus() {}

	bool intersect(const Ray& ray, glm::vec3& point, glm::vec3& normal) {
		return (glm::intersectRaySphere(ray.p, ray.d, position, 1.0, point, normal));
	}

	float sdf(const glm::vec3& p) {
		// glm::vec3 p1 = glm::inverse(sca) * glm::vec4(p, 1);
		// took too long to reasonably run, unsure if efficiency can be improved.

		glm::vec2 q = glm::vec2(glm::length(glm::vec2(position.x - p.x, position.z - p.z)) - t.x, position.y - p.y);
		return glm::length(q) - t.y;
	}

	glm::vec3 getNormal(const glm::vec3& p) { 
		return glm::normalize(p - position); 
	}

	void draw() {
		ofDrawSphere(position, 1.0);
	}

	glm::vec2 t = glm::vec2(1.0, 0.25);

	glm::mat4 trans;
	glm::mat4 rot1;
	glm::mat4 rot2;
	glm::mat4 sca;
};

/*
* Primitive Box Class
*/
class Box : public SceneObject {
public:
	Box(glm::vec3 p, float base, ofColor diffuse = ofColor::lightGray) { position = p; b = base; diffuseColor = diffuse; }
	Box() {}

	bool intersect(const Ray& ray, glm::vec3& point, glm::vec3& normal) {
		return (glm::intersectRaySphere(ray.p, ray.d, position, b * 2, point, normal));
	}

	float sdf(const glm::vec3& p) {
		glm::vec3 q = glm::abs(position - p) - b;
		return glm::length(max(q, glm::vec3(0,0,0))) + min(max(q.x, max(q.y, q.z)), 0.0f);
	}

	glm::vec3 getNormal(const glm::vec3& p) {
		return glm::normalize(p - position);
	}

	void draw() {
		ofDrawBox(position, b);
	}

	float b = 1.0;
};

/*
* Basic 3D Mandelbulb fractal, looks like a pinecone
* Source: https://iquilezles.org/articles/mandelbulb/
*/
class Mandelbulb : public SceneObject {
public:
	Mandelbulb(glm::vec3 p, int pwr = 8, int iterations = 10, float sca = 2, ofColor diffuse = ofColor::lightBlue) {
		position = p;
		power = pwr;
		iter = iterations;
		scale = sca;
		diffuseColor = diffuse;
	}
	Mandelbulb() {}

	bool intersect(const Ray& ray, glm::vec3& point, glm::vec3& normal) {
		return (glm::intersectRaySphere(ray.p, ray.d, position, radius, point, normal));
	}

	/*
	* Distance function for basic mandelbulb shape which looks like a pinecone.
	* Utilizes the basic z^p + c algorithm
	*/
	float sdf(const glm::vec3& p) {
		// rotation matrix for object
		glm::mat3x3 rot = glm::mat3x3(
			glm::vec3(cos(0.5), -sin(0.5), 0.0),
			glm::vec3(sin(0.5), cos(0.5), 0.0),
			glm::vec3(0.0, 0.0, 1.0)
		);

		glm::vec3 offset = rot * (position - p);
		glm::vec3 z = offset;
		float m = dot(z, z);

		float dz = 1.0;

		for (int i = 0; i < iter; i++)
		{
			dz = power * pow(m, 3.5) * dz + 1.0;

			float r = length(z);
			float b = power * acos(z.y / r);
			float a = power * atan(z.x / z.z);

			// z^p + c
			z = offset + pow(r, power) * glm::vec3(sin(b) * sin(a), cos(b), sin(b) * cos(a));       


			m = dot(z, z);
			if (m > bp)
				break;
		}

		// Hubbard-Douady potential
		return 0.25 * log(m) * sqrt(m) / dz;
	}

	// unused for rendering, getNormalRM is used instead.
	glm::vec3 getNormal(const glm::vec3& p) { return glm::normalize(p - position); }

	/*
	* Method to draw temporary sphere in the viewport.
	* This will not be shown in the Raymarched render.
	*/
	void draw() {
		//	ofDrawSphere(position, radius); 
		spherePrim.setRadius(radius);
		spherePrim.setPosition(position);
		spherePrim.draw();
	}

	float radius = 1.0;
	ofSpherePrimitive spherePrim;
	int power = 8;
	int iter = 10;
	float scale = 2.0;
	float bp = 256.0;
};

/*
* A 3D Representation of a basic Julia Set, the shape is a recursive tetrahedron
*/
class Julia : public SceneObject {
public:
	Julia(glm::vec3 p, int iterations = 10, float sca = 2, ofColor diffuse = ofColor::yellow) { 
		position = p; 
		iter = iterations;
		scale = sca;
		diffuseColor = diffuse; }
	Julia() {}
	bool intersect(const Ray& ray, glm::vec3& point, glm::vec3& normal) {
		return (glm::intersectRaySphere(ray.p, ray.d, position, radius, point, normal));
	}

	/*
	* Distance function for the 3D Julia Set, which is a recursive tetrahedron.
	*/
	float sdf(const glm::vec3& p) {
		// rotation matrix for object
		glm::mat3x3 rot = glm::mat3x3(
			glm::vec3(cos(0.5), -sin(0.5), 0.0),
			glm::vec3(sin(0.5), cos(0.5), 0.0),
			glm::vec3(0.0, 0.0, 1.0)
		);

		glm::vec3 a1 = rot * glm::vec3( 1.0,  1.0,  1.0);
		glm::vec3 a2 = rot * glm::vec3(-1.0,  1.0, -1.0);
		glm::vec3 a3 = rot * glm::vec3( 1.0, -1.0, -1.0);
		glm::vec3 a4 = rot * glm::vec3(-1.0, -1.0,  1.0);

		glm::vec3 z = position - p;
		glm::vec3 c;
		int n = 0;
		float dist, d;
		while (n < iter) {
			c = a1; dist = length(z - a1);
			d = length(z - a2); 
			if (d < dist) 
			{ 
				c = a2; 
				dist = d; 
			}
			d = length(z - a3); 
			if (d < dist) {
				c = a3; 
				dist = d; 
			}
			d = length(z - a4); 
			if (d < dist) { 
				c = a4; 
				dist = d; 
			}

			//z^p - c
			z = scale * z - c * (scale - 1.0);
			n++;
		}

		return length(z) * pow(scale, float(-n)) * 0.5;
	}

	// unused for rendering, getNormalRM is used instead.
	glm::vec3 getNormal(const glm::vec3& p) { return glm::normalize(p - position); }

	/*
	* Method to draw temporary sphere in the viewport.
	* This will not be shown in the Raymarched render.
	*/
	void draw() {
		//	ofDrawSphere(position, radius); 
		spherePrim.setRadius(radius);
		spherePrim.setPosition(position);
		spherePrim.draw();
	}

	float radius = 1.0;
	ofSpherePrimitive spherePrim;
	int iter = 10;
	float scale = 2.0;
};

//  General purpose plane 
//
class Plane: public SceneObject {
public:
	Plane(glm::vec3 p, glm::vec3 n, ofColor diffuse = ofColor::green, float w = 20, float h = 20 ) {
		position = p;
		normal = n;
		width = w;
		height = h;
		diffuseColor = diffuse;
		if (normal == glm::vec3(0, 1, 0))
			plane.rotateDeg(-90, 1, 0, 0);
		else if (normal == glm::vec3(0, -1, 0))
			plane.rotateDeg(90, 1, 0, 0);
		else if (normal == glm::vec3(1, 0, 0))
			plane.rotateDeg(90, 0, 1, 0);
		else if (normal == glm::vec3(-1, 0, 0))
			plane.rotateDeg(-90, 0, 1, 0);
	}
	Plane() { 
		normal = glm::vec3(0, 1, 0);
		plane.rotateDeg(90, 1, 0, 0);
		isSelectable = false;
	}

	bool intersect(const Ray &ray, glm::vec3 & point, glm::vec3 & normal);

	float sdf(const glm::vec3 & p);
	// return glm::dot(p, getNormal(p)) + height;

	glm::vec3 getNormal(const glm::vec3 &p) { return this->normal; }
	void draw() {
		plane.setPosition(position);
		plane.setWidth(width);
		plane.setHeight(height);
		plane.setResolution(4, 4);
	//	plane.drawWireframe();
		plane.draw();
	}
	ofPlanePrimitive plane;
	glm::vec3 normal;
	
	float width = 20;
	float height = 20;

};

// view plane for render camera
// 
class  ViewPlane: public Plane {
public:
	ViewPlane(glm::vec2 p0, glm::vec2 p1) { min = p0; max = p1; }

	ViewPlane() {                         // create reasonable defaults (6x4 aspect)
		min = glm::vec2(-3, -2);
		max = glm::vec2(3, 2);
		position = glm::vec3(0, 0, 5);
		normal = glm::vec3(0, 0, 1);      // viewplane currently limited to Z axis orientation
	}

	void setSize(glm::vec2 min, glm::vec2 max) { this->min = min; this->max = max; }
	float getAspect() { return width() / height(); }

	glm::vec3 toWorld(float u, float v);   //   (u, v) --> (x, y, z) [ world space ]

	void draw() {
		ofDrawRectangle(glm::vec3(min.x, min.y, position.z), width(), height());
	}

	
	float width() {
		return (max.x - min.x);
	}
	float height() {
		return (max.y - min.y); 
	}

	// some convenience methods for returning the corners
	//
	glm::vec2 topLeft() { return glm::vec2(min.x, max.y); }
	glm::vec2 topRight() { return max; }
	glm::vec2 bottomLeft() { return min;  }
	glm::vec2 bottomRight() { return glm::vec2(max.x, min.y); }

	//  To define an infinite plane, we just need a point and normal.
	//  The ViewPlane is a finite plane so we need to define the boundaries.
	//  We will define this in terms of min, max  in 2D.  
	//  (in local 2D space of the plane)
	//  ultimately, will want to locate the ViewPlane with RenderCam anywhere
	//  in the scene, so it is easier to define the View rectangle in a local'
	//  coordinate system.
	//
	glm::vec2 min, max;

	
};



//  Light Base Class
//

class Light: public SceneObject {
public:
	float intensity = 1.0;

	// intersection routine for UI (not ray tracing) - used to move lights with mouse
	//
	bool mouseIntersect(const Ray &ray, glm::vec3 &point, glm::vec3 &normal) {
		return (glm::intersectRaySphere(ray.p, ray.d, position, .2, point, normal));
	}

	ofLight vpLight;             // for lighting in viewport
};

//  Point Light
//
class PointLight: public Light {
public:
	PointLight(const glm::vec3 &p, float intensity = 1.0, ofColor diffuse = ofColor::white) {
		position = p;
		diffuseColor = diffuse;
		this->intensity = intensity;
		vpLight.setPosition(position);
	}
	void draw() {
		ofSetColor(diffuseColor);
		ofDrawSphere(position, .2);
	}
	
};

class AmbientLight : public Light {
public:
	AmbientLight(ofColor diffuse = ofColor::white) { diffuseColor = diffuse; intensity = .05;  }  // ambient light has no position
	void draw() {
		ofSetColor(diffuseColor);
		ofDrawSphere(position, .2);
	}
};


//  render camera  - currently must be z axis aligned (we will improve this in project 4)
//
class RenderCam: public SceneObject {
public:
	RenderCam() {
		position = glm::vec3(0, 0, 10);
		aim = glm::vec3(0, 0, -1);
	}
	Ray getRay(float u, float v);
	void draw() { ofDrawBox(position, 1.0); };
	void drawFrustum();

	glm::vec3 aim;
	ViewPlane view;          // The camera viewplane, this is the view that we will render 
};

  

class ofApp : public ofBaseApp{
	public:
		void setup();
		void update();
		void draw();

		void keyPressed(int key);
		void keyReleased(int key);
		void mouseMoved(int x, int y );
		void mouseDragged(int x, int y, int button);
		void mousePressed(int x, int y, int button);
		void mouseReleased(int x, int y, int button);
		void mouseEntered(int x, int y);
		void mouseExited(int x, int y);
		void windowResized(int w, int h);
		void dragEvent(ofDragInfo dragInfo);
		void gotMessage(ofMessage msg);
		void rayTrace();
		void rayMarchRender();
		bool rayMarch(const Ray &r, glm::vec3 & point, int & obj);
		void drawGrid();
		void drawAxis(glm::vec3 position);
		bool mouseToWorld(int x, int y, glm::vec3 &point);
		void addLight(Light *light) {
			lights.push_back(light);
		}

		ofColor lambert(const glm::vec3 &p, const glm::vec3 &n, const ofColor diffuse);
		ofColor phong(const glm::vec3 &p, const glm::vec3 &n, const ofColor diffuse, const ofColor specular, float power);
	
		// Ray Marching functions
		//
		float sceneSDF(glm::vec3 p, int& obj);
		float sceneSDF(glm::vec3 p);
		glm::vec3 getNormalRM(const glm::vec3& p);
		ofColor phongRM(const glm::vec3 &p, const glm::vec3 &n, const ofColor diffuse, const ofColor specular, float power);

		bool bHide = true;
		bool bShowImage = false;

		ofEasyCam  mainCam;
		ofCamera sideCam;
		ofCamera previewCam;
		ofCamera  *theCam;    // set to current camera either mainCam or sideCam

		// set up one render camera to render image throughn
		//
		RenderCam renderCam;
		ofImage image;
		ofImage previewImage;

		// scene components
		//
		vector<SceneObject *> scene;
		vector<Light *> lights;
		AmbientLight ambientLight;
		ofSpherePrimitive selectionMarker;

		SceneObject *selectedObj = NULL;
		glm::vec3 selectedPoint;

		int imageWidth =  1200;
		int imageHeight = 800;

		bool bDrag = false;
		bool bRayMarched = false;
		glm::vec3 lastPoint;

};