/* mouse2.c 

Laura Toma

What it does:  The user can enter a polygon by clicking on the mouse. 

*/

#include "geom.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif


#include <vector> 

using namespace std; 

GLfloat red[3] = {1.0, 0.0, 0.0};
GLfloat green[3] = {0.0, 1.0, 0.0};
GLfloat blue[3] = {0.0, 0.0, 1.0};
GLfloat black[3] = {0.0, 0.0, 0.0};
GLfloat white[3] = {1.0, 1.0, 1.0};
GLfloat gray[3] = {0.5, 0.5, 0.5};
GLfloat yellow[3] = {1.0, 1.0, 0.0};
GLfloat magenta[3] = {1.0, 0.0, 1.0};
GLfloat cyan[3] = {0.0, 1.0, 1.0};

GLint fillmode = 0;

/* forward declarations of functions */
void display(void);
void keypress(unsigned char key, int x, int y);
void mousepress(int button, int state, int x, int y);
void timerfunc(); 

void initialize_polygon(); 
void print_polygon(vector<point2D> poly);

/* our coordinate system is (0,0) x (WINDOWSIZE,WINDOWSIZE) where the
   origin is the lower left corner */

/* global variables */
const int WINDOWSIZE = 750; 

vector<eventPoint> events;

vector<point2D> visiblePolygonPoints;

vector<lineSegment2D> lines;

vector< vector<point2D> > triangulate;

//the current polygon 
vector<point2D> poly;

point2D guardPoint;

double currentDirectionX = 1;
double currentDirectionY = 1;

vector<point2D> polygonLines;

bool moveGuard = false;

//coordinates of last mouse click
double mouse_x=-10, mouse_y=-10; 
//initialized to a point outside the window

//when this is 1, then clicking the mouse results in those points being stored in poly
int poly_init_mode = 0; 

void draw_circle(double x, double y){
  glColor3fv(blue);   
  int precision = 100;
  double r = 4; 
  double theta = 0;
  glBegin(GL_POLYGON);
  for(int i = 0; i < precision; i++){
    theta = i * 2 * M_PI/precision;
    glVertex2f(x + r*cos(theta), y + r*sin(theta));
  }
  glEnd();
}

/* 
Usage

void glutMouseFunc(void (*func)(int button, int state, int x, int y));

Description

glutMouseFunc sets the mouse callback for the current window. When a
user presses and releases mouse buttons in the window, each press and
each release generates a mouse callback. The button parameter is one
of GLUT_LEFT_BUTTON, GLUT_MIDDLE_BUTTON, or GLUT_RIGHT_BUTTON. For
systems with only two mouse buttons, it may not be possible to
generate GLUT_MIDDLE_BUTTON callback. For systems with a single mouse
button, it may be possible to generate only a GLUT_LEFT_BUTTON
callback. The state parameter is either GLUT_UP or GLUT_DOWN
indicating whether the callback was due to a release or press
respectively. The x and y callback parameters indicate the window
relative coordinates when the mouse button state changed. If a
GLUT_DOWN callback for a specific button is triggered, the program can
assume a GLUT_UP callback for the same button will be generated
(assuming the window still has a mouse callback registered) when the
mouse button is released even if the mouse has moved outside the
window.
*/

void mousepress(int button, int state, int x, int y) {


  if (state == GLUT_DOWN) {

    mouse_x = x;
    mouse_y = y;
    //(x,y) are in wndow coordinates, where the origin is in the upper
    //left corner; our reference system has the origin in lower left
    //corner, this means we have to reflect y
    mouse_y = WINDOWSIZE - mouse_y; 

    printf("mouse click at (x=%d, y=%d)\n", (int)mouse_x, (int)mouse_y);

    if (poly_init_mode == 1) {

      triangulate.clear();
      visiblePolygonPoints.clear();

      point2D p = {mouse_x, mouse_y}; 
      poly.push_back(p);
    }

    else {
      guardPoint.x = mouse_x;
      guardPoint.y = mouse_y;

      if (!isSimple(poly)) {
        cout << "Polygon is not simple" << endl;
        triangulate.clear();
      }

      else if (!isInPolygon(poly, guardPoint)) {
        cout << "Point is not in the polygon" << endl;
        triangulate.clear();
      }

      else {
        print_polygon(poly);
        triangulate = triangulatedPolygon(poly, guardPoint);
      }

    }

  }
  
  glutPostRedisplay();
}

/* ****************************** */
/* initialize  polygon stored in global variable poly  */
void initialize_polygon() {
  
  //clear the vector, in case something was there 
  poly.clear(); 

  int n = 10; //size of polygon 
  double rad = 100; 
  double step = 2 * M_PI / n;
  point2D p; 
  for (int i=0; i<n; i++) {

    p.x = WINDOWSIZE/2 + rad * cos (i * step); 
    p.y = WINDOWSIZE/2 + rad * sin (i * step); 

    //insert the segment in the array of segments 
    poly.push_back (p); 
  } //for i
}

/* ************************************************** */
void print_polygon(vector<point2D> poly) {

  if (poly.size() == 0 ) {
    return;
  }

  for (unsigned int i=0; i<poly.size()-1; i++) {
    printf("edge %d: [(%d,%d), (%d,%d)]\n",
	   i, poly[i].x, poly[i].y, poly[i+1].x, poly[i+1].y);
  }
  //print last edge from last point to first point 
  int last = poly.size()-1; 
    printf("edge %d: [(%d,%d), (%d,%d)]\n",
	   last, poly[last].x, poly[last].y, poly[0].x, poly[0].y);

}

/* ****************************** */
int main(int argc, char** argv) {

  initialize_polygon();
  print_polygon(poly);


  /* initialize GLUT  */
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
  glutInitWindowSize(WINDOWSIZE, WINDOWSIZE);
  glutInitWindowPosition(100,100);
  glutCreateWindow(argv[0]);

  /* register callback functions */
  glutDisplayFunc(display); 
  glutKeyboardFunc(keypress);
  glutMouseFunc(mousepress); 
  glutIdleFunc(timerfunc); 

  /* init GL */
  /* set background color black*/
  glClearColor(0, 0, 0, 0);   
  /* here we can enable depth testing and double buffering and so
     on */
  
  /* give control to event handler */
  glutMainLoop();
  return 0;
}

/* ****************************** */
/* draw the polygon */
void draw_polygon(vector<point2D> poly){

  if (poly.size() == 0) return; 

  glColor3fv(yellow); 

  int i;
  for (i=0; i<poly.size()-1; i++) {

    glBegin(GL_LINES);
    glVertex2f(poly[i].x, poly[i].y); 
    glVertex2f(poly[i+1].x, poly[i+1].y);
    glEnd();
  }
  //render last segment between last point and forst point 
  int last=poly.size()-1; 
    glBegin(GL_LINES);
    glVertex2f(poly[last].x, poly[last].y); 
    glVertex2f(poly[0].x, poly[0].y);
    glEnd();
}

void render_visible_polygon() {

  glColor3fv(red);

  for (int i = 0; i < triangulate.size(); i++) {

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glBegin(GL_TRIANGLES);
    glVertex2f(triangulate[i][0].x, triangulate[i][0].y);
    glVertex2f(triangulate[i][1].x, triangulate[i][1].y);
    glVertex2f(triangulate[i][2].x, triangulate[i][2].y);
    glEnd();

    }

}

/* ****************************** */
void display(void) {

  glClear(GL_COLOR_BUFFER_BIT);
  glMatrixMode(GL_MODELVIEW); 
  glLoadIdentity(); //clear the matrix


  /* The default GL window is [-1,1]x[-1,1] with the origin in the
     center. 
     
     Our system of coordinates (in which we generate our points) is
     (0,0) to (WINSIZE,WINSIZE), with the origin in the lower left
     corner.
     
     We need to map the points to [-1,1] x [-1,1]  
     
     Assume we are the local coordinate system. 
     
     First we scale down to [0,2] x [0,2] */ 
  glScalef(2.0/WINDOWSIZE, 2.0/WINDOWSIZE, 1.0);  
   /* Then we translate so the local origin goes in the middle of teh
     window to (-WINDOWSIZE/2, -WINDOWSIZE/2) */
  glTranslatef(-WINDOWSIZE/2, -WINDOWSIZE/2, 0); 

  //now we draw in our local coordinate system (0,0) to
  //(WINSIZE,WINSIZE), with the origin in the lower left corner.

  render_visible_polygon();

  draw_polygon(poly);
  
  //draw a circle where the mouse was last clicked. Note that this
  //point is stored as a global variable and is modified by the mouse
  //handler function
  draw_circle(mouse_x, mouse_y);

  /* execute the drawing commands */
  glFlush();
}

/* ****************************** */
void keypress(unsigned char key, int x, int y) {
  switch(key) {
  case 'q':
    exit(0);
    break;

    //expected behaviour: press 's', then click on the points you
    //want, and press 'e' when you're done. the points will be saved
    //in global 'poly'
    
  case 's': 
    //start re-initializeing polygon 
    poly.clear();
    events.clear();
    visiblePolygonPoints.clear();
    lines.clear();
    triangulate.clear();
    mouse_x = mouse_y=0; 
    poly_init_mode = 1; 
    glutPostRedisplay();
    break; 
    
  case 'e':
    poly_init_mode = 0; 
    glutPostRedisplay();
    break;  

  case 'm':
    moveGuard = true;
    glutPostRedisplay();
    break;


  case 'n':
    moveGuard = false;
    glutPostRedisplay();
    break;

  }

}

/* Handler for window re-size event. Called back when the window first appears and
   whenever the window is re-sized with its new width and height */
void reshape(GLsizei width, GLsizei height) {  // GLsizei for non-negative integer
     
   // Set the viewport to cover the new window
   glViewport(0, 0, width, height);
 
   glMatrixMode(GL_PROJECTION);  // To operate on the Projection matrix
   glLoadIdentity();             // Reset
   gluOrtho2D(0.0, (GLdouble) width, 0.0, (GLdouble) height); 
}

bool currDirectionValid() {

  point2D testPoint;
  testPoint.x = mouse_x += currentDirectionX;
  testPoint.y = mouse_y += currentDirectionY;

  if (isInPolygon(poly, testPoint)) {
    return true;
  }

  triangulate.clear();
  return false;

}

int randomDirection() {
  return (rand() % 3) - 1;
}

void determineNewDirection() {

  while (!currDirectionValid()) {
    triangulate.clear();
    currentDirectionX = randomDirection();
    currentDirectionY = randomDirection();
  }

}

void timerfunc() {

  if (moveGuard == true) {

    if (currDirectionValid()) {
      guardPoint.x += currentDirectionX;
      guardPoint.y += currentDirectionY;
    }

    else {
      determineNewDirection();
    }

    triangulate = triangulatedPolygon(poly, guardPoint);
    glutPostRedisplay();

  }


}
