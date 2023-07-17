#include "ofApp.h"

// Intersect Ray with Plane  (wrapper on glm::intersect*)
//
bool Plane::intersect(const Ray& ray, glm::vec3& point, glm::vec3& normalAtIntersect) {
	float dist;
	bool insidePlane = false;
	bool hit = glm::intersectRayPlane(ray.p, ray.d, position, this->normal, dist);
	if (hit) {
		Ray r = ray;
		point = r.evalPoint(dist);
		normalAtIntersect = this->normal;
		glm::vec2 xrange = glm::vec2(position.x - width / 2, position.x + width / 2);
		glm::vec2 yrange = glm::vec2(position.y - width / 2, position.y + width / 2);
		glm::vec2 zrange = glm::vec2(position.z - height / 2, position.z + height / 2);

		// horizontal 
		//
		if (normal == glm::vec3(0, 1, 0) || normal == glm::vec3(0, -1, 0)) {
			if (point.x < xrange[1] && point.x > xrange[0] && point.z < zrange[1] && point.z > zrange[0]) {
				insidePlane = true;
			}
		}
		// front or back
		//
		else if (normal == glm::vec3(0, 0, 1) || normal == glm::vec3(0, 0, -1)) {
			if (point.x < xrange[1] && point.x > xrange[0] && point.y < yrange[1] && point.y > yrange[0]) {
				insidePlane = true;
			}
		}
		// left or right
		//
		else if (normal == glm::vec3(1, 0, 0) || normal == glm::vec3(-1, 0, 0)) {
			if (point.y < yrange[1] && point.y > yrange[0] && point.z < zrange[1] && point.z > zrange[0]) {
				insidePlane = true;
			}
		}
	}
	return insidePlane;
}

float Plane::sdf(const glm::vec3 & p) {
	// for the moment, lets just consider ground plane
	//
    //	float freq = .07;
    //	float amp = 3.0;
	if (normal == glm::vec3(0, 1, 0)) {
		return p.y - position.y;
	}
	else if (normal == glm::vec3(0, 0, 1)) {
		return p.z - position.z;
	}
	else {
		return 0.0;
	}
}

// Convert (u, v) to (x, y, z) 
// We assume u,v is in [0, 1]
//
glm::vec3 ViewPlane::toWorld(float u, float v) {
	float w = width();
	float h = height();
	return (glm::vec3((u * w) + min.x, (v * h) + min.y, position.z));
}

// Get a ray from the current camera position to the (u, v) position on
// the ViewPlane
//
Ray RenderCam::getRay(float u, float v) {
	glm::vec3 pointOnPlane = view.toWorld(u, v);
	return(Ray(position, glm::normalize(pointOnPlane - position)));
}

// This could be drawn a lot simpler but I wanted to use the getRay call
// to test it at the corners.
// 
void RenderCam::drawFrustum() {
	view.draw();
	Ray r1 = getRay(0, 0);
	Ray r2 = getRay(0, 1);
	Ray r3 = getRay(1, 1);
	Ray r4 = getRay(1, 0);
	float dist = glm::length((view.toWorld(0, 0) - position));
	r1.draw(dist);
	r2.draw(dist);
	r3.draw(dist);
	r4.draw(dist);
}

//--------------------------------------------------------------
void ofApp::setup() {
	ofSetBackgroundColor(ofColor::black);
	ofEnableDepthTest();
	mainCam.setDistance(15);
	mainCam.setNearClip(.1);
	sideCam.setPosition(40, 0, 0);
	sideCam.lookAt(glm::vec3(0, 0, 0));
	previewCam.setPosition(renderCam.position);
	previewCam.setNearClip(.1);
	previewCam.lookAt(glm::vec3(0, 0, 0));
	theCam = &mainCam;

	
	// create a scene	
	// 3D Primitives
	scene.push_back(new Torus(glm::vec3(-3, -1, 0), ofColor::pink));
	scene.push_back(new Box(glm::vec3(-3, -1, 0), 0.5, ofColor::blue)); //
	//scene.push_back(new Sphere(glm::vec3(-2, 0, 1), 0.5, ofColor::orange));
	//scene.push_back(new Sphere(glm::vec3(-1, 1, 0.5), 0.5, ofColor::yellow));
    //scene.push_back(new Sphere(glm::vec3(0, 2, 0), 0.5, ofColor::green));
	//scene.push_back(new Sphere(glm::vec3(1, 1, 0.5), 0.5, ofColor::blue));
	//scene.push_back(new Sphere(glm::vec3(2, 0, 1), 0.5, ofColor::indigo));
	//scene.push_back(new Box(glm::vec3(3, -1, 0), 0.5, ofColor::purple)); //
	//scene.push_back(new Torus(glm::vec3(3, -1, 0), ofColor::pink));

	// 3D Fractal Objects
	//scene.push_back(new Mandelbulb(glm::vec3(-3, 1, 0), 8, 8, 2.0));
	scene.push_back(new Julia(glm::vec3(3, 1, 0), 8, 2.0));

	// Ground plane
	scene.push_back(new Plane(glm::vec3(0, -2, 0), glm::vec3(0, 1, 0), ofColor::brown));

	// Scene Lights
	addLight(new PointLight(glm::vec3(0, 10, 0), 100));
	//addLight(new PointLight(glm::vec3(10, 0, 0), 100)); // 50

	// output image to render to
	image.allocate(imageWidth, imageHeight, ofImageType::OF_IMAGE_COLOR);
	previewImage.allocate(imageWidth, imageHeight, ofImageType::OF_IMAGE_COLOR);
	
}

// main ray trace loop (for Assignment Part 1) - Kevin M. Smith  10/4/2018
// not implemeneted
void ofApp::rayTrace() {

	// iterate in "pixel" space starting from lower-left corner
	//
	for (int j = 0; j < imageHeight; j++) {
		for (int i = 0; i < imageWidth; i++) {

		}
	}
}

/* Main ray marching loop - Calvin Quach 3/12/2023
*  For each pixel, march the ray to check for a surface hit.
*/
void ofApp::rayMarchRender() {

	// iterate in "pixel" space starting from lower-left corner
	//
	for (int j = 0; j < imageHeight; j++) {
		for (int i = 0; i < imageWidth; i++) {

			// from slides from lecture; compute a normalized (u, v) in [0  1.0]
			// from the pixel coordinates. we get back a (u, v) that is centered
			//  in pixel(i, j). "v" is up.  The "float" casting shoud not be 
			// necessary for newer C++ compilers, but I'm just being explicit here.
			//
			float u = (float(i) + .5) / float(imageWidth);
			float v = (float(j) + .5) / float(imageHeight);

			// we have a helper function that gives world space coordinates (x, y, z)
			// from (u, v)
			//
			glm::vec3 point = renderCam.view.toWorld(u, v);

		    //  March along ray to see if we have a hit
			//
			Ray r = renderCam.getRay(u, v);
			int obj = -1;
			bool hit = rayMarch(r, point, obj);
			
			// Set the color of the pixel to the nearest object's pixel
			// if we didn't hit anything, set it the bg color (Black for now).
			//
			// This is a little tricky. Some image file formats and API's
			// assume that the pixel (0, 0) is in the upper left and this is the 
			// case with OpenFrameworks jpeg images.  So we need to make sure
			// we adjust for that in the "j" direction, otherwise our image
			// will come out "flipped".
			//
			if (!hit) {
				image.setColor(i, imageHeight - j - 1, ofColor::black);
			}
			else {
				// Use ray marching version of phong shader
				// 
				ofColor color = phongRM(point, getNormalRM(point), scene[obj]->diffuseColor, scene[obj]->specularColor, 10.0);
				image.setColor(i, imageHeight - j - 1, color);
			}
		}
	}
	image.save("out-rm.jpg");
}


/* Main ray march algorithm
*  The marching ray has hit a surface if its distance 
*  to point p is within the distance threshold.
*  Calvin Quach 3/12/2023
*/
bool ofApp::rayMarch(const Ray& r, glm::vec3& point, int& obj) {
	bool hit = false;
	glm::vec3 p = r.p;
	for (int i = 0; i < MAX_RAY_STEPS; i++)
	{
		float dist = sceneSDF(p, obj);
		if (dist < DIST_THRESHOLD)
		{
			hit = true;
			point = p;
			break;
		}
		else if (dist > MAX_DISTANCE) {
			break;
		}
		else {
			p = p + r.d * dist; // to move along ray
		}
	}

	return hit;
}

/* Method to calculate the closest object along with its distnace at the respective point p.
*  Calvin Quach 3/12/2023
*/
float ofApp::sceneSDF(glm::vec3 p, int& obj)
{
	float closestDis = std::numeric_limits<float>::infinity();
	for (int i = 0; i < scene.size(); i++)
	{
		float d = scene[i]->sdf(p);
		if (d < closestDis)
		{
			closestDis = d;
			obj = i;
		}
	}
	return closestDis;
}

/* Method to calculate the closest distnace at the respective point p.
*  Calvin Quach 3/12/2023
*/
float ofApp::sceneSDF(glm::vec3 p)
{
	float closestDis = std::numeric_limits<float>::infinity();
	for (int i = 0; i < scene.size(); i++)
	{
		float d = scene[i]->sdf(p);
		if (d < closestDis)
		{
			closestDis = d;
		}
	}
	return closestDis;
}

/* Method to get the gradient normal of the respective point 
*  Calvin Quach 3/12/2023
*/
glm::vec3 ofApp::getNormalRM(const glm::vec3& p) {
	float eps = 0.01; // .0001 or .01
	float dp = sceneSDF(p);
	glm::vec3 n = glm::vec3(
		dp - sceneSDF(glm::vec3(p.x - eps, p.y, p.z)),
		dp - sceneSDF(glm::vec3(p.x, p.y - eps, p.z)),
		dp - sceneSDF(glm::vec3(p.x, p.y, p.z - eps))
	);

	return glm::normalize(n);
}

//  Lambert shader
//
ofColor ofApp::lambert(const glm::vec3 &p, const glm::vec3 &norm, const ofColor diffuse) {
	ofColor color;
	return color;
}

//  Phong shader (includes lambert + phone/blinn component)
//
ofColor ofApp::phong(const glm::vec3 &p, const glm::vec3 &norm, const ofColor diffuse, const ofColor specular, float power) {
	ofColor color;
	return color;
}

/* Phong shader helper method (includes lambert + phone/blinn component)
*  This version is adapated for Ray Marching.
*  Calvin Quach 3/12/2023
*/
ofColor ofApp::phongRM(const glm::vec3 &p, const glm::vec3 &norm, const ofColor diffuse, const ofColor specular, float power) {
	ofColor ambient = diffuse * 0.1;
	for (int i = 0; i < lights.size(); i++)
	{
		glm::vec3 n = norm;
		glm::vec3 v = glm::normalize(renderCam.position - p);
		glm::vec3 l = glm::normalize(lights[i]->position - p);
		glm::vec3 h = glm::normalize(v + l);
		float r = glm::distance(p, lights[i]->position);

		int obj = -1;
		glm::vec3 po = p;
		po.y += 0.1;

		if (!rayMarch(Ray(po, l), po, obj))
		{
			ofColor diff = diffuse;
			ofColor spec = specular;

			ambient += diff * (lights[i]->intensity / glm::pow(r, 2)) * glm::max(0.0f, glm::dot(n, l));
			ambient += spec * (lights[i]->intensity / glm::pow(r, 2)) * glm::max(0.0f, glm::pow(glm::dot(n, h), power));
		}
	}

	return ambient;
}

//--------------------------------------------------------------
void ofApp::drawGrid() {

	float u = 0;
	float v = 0;
	float pixelWidth = 1.0 / imageWidth;
	float pixelHeight = 1.0 / imageHeight;
	for (int x = 0; x < imageWidth; x++) {
		glm::vec3 p1 = renderCam.view.toWorld(u, 0);
		glm::vec3 p2 = renderCam.view.toWorld(u, 1);
		ofDrawLine(p1, p2);
		u += pixelWidth;
	}
	for (int y = 0; y < imageHeight; y++) {
		glm::vec3 p1 = renderCam.view.toWorld(0, v);
		glm::vec3 p2 = renderCam.view.toWorld(1, v);
		ofDrawLine(p1, p2);
		v += pixelHeight;
	}
}

 
//--------------------------------------------------------------
void ofApp::update(){

}

//--------------------------------------------------------------
void ofApp::draw(){

	//ofEnableLighting();

	// enable viewport lighting
	//
//	for (int i = 0; i < lights.size(); i++) {
//		lights[i]->vpLight.enable();
//	}


	theCam->begin();
	
	ofSetColor(ofColor::green);
	ofNoFill();

	drawAxis(glm::vec3(0, 0, 0));;

	//  draw objects in scene
	//
	for (int i = 0; i < scene.size(); i++) {
		if (scene[i] == selectedObj)
			ofSetColor(ofColor::grey);
		else ofSetColor(scene[i]->diffuseColor);
		scene[i]->draw();
	}

	if (selectedObj) {
		selectionMarker.setPosition(selectedPoint);
		selectionMarker.setRadius(.1);
		ofSetColor(ofColor::yellow);
		selectionMarker.draw();
	}

	ofPushMatrix();
	glm::vec3 eyePos = glm::vec3(0.001, 5.0, 0);
	glm::vec3 aimPos = glm::vec3(0, 0, 0);
	glm::vec3 upVector = glm::vec3(0, 1, 0);
	glm::mat4 m = glm::lookAt(eyePos, aimPos, upVector);
	ofMultMatrix(glm::inverse(m));
	ofRotate(-90, 1, 0, 0);
	ofSetColor(ofColor::lightGray);
	ofDrawCone(.5, 1);
	ofPopMatrix();

	ofDisableLighting();
	ofSetColor(ofColor::lightSkyBlue);
	renderCam.drawFrustum();
	ofSetColor(ofColor::blue);
	renderCam.draw();
//	drawGrid();

	theCam->end();

	if (bShowImage) {
		ofSetColor(255, 255, 255);
		previewImage.draw((ofGetWindowWidth() - image.getWidth()) / 2, (ofGetWindowHeight() - image.getHeight()) / 2);
	}
}

// 
// Draw an XYZ axis in RGB at world (0,0,0) for reference.
//
void ofApp::drawAxis(glm::vec3 position) {
	ofPushMatrix();
	ofTranslate(position);

	ofSetLineWidth(1.0);

	// X Axis
	ofSetColor(ofColor(255, 0, 0));
	ofDrawLine(ofPoint(0, 0, 0), ofPoint(1, 0, 0));


	// Y Axis
	ofSetColor(ofColor(0, 255, 0));
	ofDrawLine(ofPoint(0, 0, 0), ofPoint(0, 1, 0));

	// Z Axis
	ofSetColor(ofColor(0, 0, 255));
	ofDrawLine(ofPoint(0, 0, 0), ofPoint(0, 0, 1));

	ofPopMatrix();
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key){

}

//--------------------------------------------------------------
void ofApp::keyReleased(int key) {
	switch (key) {
	case 'C':
	case 'c':
		if (mainCam.getMouseInputEnabled()) mainCam.disableMouseInput();
		else mainCam.enableMouseInput();
		break;
	case 'F':
	case 'b':
		break;
	case 'f':
		ofToggleFullscreen();
		break;
	case 'h':
		bHide = !bHide;
		break;
	case 'i':
		bShowImage = !bShowImage;
		if (bShowImage) {
			if (bRayMarched)
				previewImage.load("out-rm.jpg");
			else
				previewImage.load("out-rt.jpg");
		}
		break;
	case 'n':
		scene.push_back(new Sphere(glm::vec3(0, 0, 0), 1.0, ofColor::violet));
		break;
	case 'r':
		bRayMarched = false;
		cout << "rendering..." << endl;
		rayTrace();
		cout << "done..." << endl;
		break;
	case 'm':
		bRayMarched = true;
		cout << "ray marching..." << endl;
		rayMarchRender();
		cout << "done..." << endl;
		break;
	case OF_KEY_F1: 
		theCam = &mainCam;
		break;
	case OF_KEY_F2:
		theCam = &sideCam;
		break;
	case OF_KEY_F3:
		theCam = &previewCam;
		break;
	}
}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y ){

}
//--------------------------------------------------------------


//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button) {

	if (mainCam.getMouseInputEnabled()) return;

	if (selectedObj && bDrag) {
		glm::vec3 point;
		mouseToWorld(x, y, point);
		selectedObj->position += (point - lastPoint);
		lastPoint = point;
	}
}

bool ofApp::mouseToWorld(int x, int y, glm::vec3 &point) {
	glm::vec3 p = theCam->screenToWorld(glm::vec3(x, y, 0));
	glm::vec3 d = p - theCam->getPosition();
	glm::vec3 dn = glm::normalize(d);

	float dist;
	if (glm::intersectRayPlane(p, d, selectedObj->position, glm::normalize(theCam->getZAxis()), dist)) {
		point = p + d * dist;
		return true;
	}
	return false;
}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button){

	if (mainCam.getMouseInputEnabled()) return;
	//
	// test if something selected
	//
	vector<SceneObject *> selected;
	selectedObj = NULL;

	glm::vec3 p = theCam->screenToWorld(glm::vec3(x, y, 0));
	glm::vec3 d = p - theCam->getPosition();
	glm::vec3 dn = glm::normalize(d);

	// check for selection of scene objects
	//
	for (int i = 0; i < scene.size(); i++) {
		
		glm::vec3 point, norm;
		
		//  We hit an object
		//
		if (scene[i]->isSelectable && scene[i]->intersect(Ray(p, dn), point, norm)) {
			selected.push_back(scene[i]);
			selectedObj = scene[i];
			selectedPoint = point;
		}
	}

	// check for selection of lights as a separate case
	//
	glm::vec3 point, norm;
	for (int n = 0; n < lights.size(); n++) {
		if (lights[n]->mouseIntersect(Ray(p, dn), point, norm)) {
			selected.push_back(lights[n]);
			selectedObj = lights[n];
		}
	}

	// if we selected more than one, pick nearest
	//
	if (selected.size() > 1) {
		float nearestDist = std::numeric_limits<float>::infinity();
		SceneObject *nearestObj = NULL;
		for (int n = 0; n < selected.size(); n++) {
			float dist = glm::length(selected[n]->position - theCam->getPosition());
			if (dist < nearestDist) {
				nearestDist = dist;
				selectedObj = selected[n];
			}	
		}
	}
	if (selected.size() == 0) {
		selectedObj = NULL;
	}
	else {
		bDrag = true;
		mouseToWorld(x, y, lastPoint);
	}
}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button){
	bDrag = false;
}


//
//--------------------------------------------------------------
void ofApp::mouseEntered(int x, int y){

}

//--------------------------------------------------------------
void ofApp::mouseExited(int x, int y){

}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h){

}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg){

}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo){ 

}

