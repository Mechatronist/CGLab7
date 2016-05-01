/* 
   Laura Toma
*/

#ifndef __geom_h
#define __geom_h

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>

#include <vector>
#include <map>
#include <utility>

using namespace std;

typedef struct _point2d {
  int x,y;
  double guardAngle;
  double distanceFromGuard; 
  bool left; // true/1 = left, false/0 = right. 
  int intersectionType; // 1 = proper, 0 = improper
  int intersectionEdge;
  int distanceFromPreviousPoint;
} point2D;

typedef struct _eventPoint {
	point2D point;
	int intersectionType; // 1 = proper, 0 = improper
	double distanceFromGuard;
  int intersectionEdge;
  int distanceFromPreviousPoint;
} eventPoint;

typedef struct _lineSegment2D {
    point2D p1, p2;
} lineSegment2D;


typedef struct _rect2D  {
    point2D origin;
    float width, height;
} rect2D;

//vector<lineSegment2D> initializePolygonEdges(vector<point2D> polygonPoints);
vector<point2D> polygonShift(vector<point2D> originalPolygonPoints, point2D guardPoint);
double computeXIntersection (int firstPointX, int firstPointY, int secondPointX, int secondPointY);
bool isInPolygon(vector<point2D> polygonPoints, point2D guardPoint);
bool collinear(point2D a, point2D b, point2D c);
vector<lineSegment2D> makePolygonEdges(vector<point2D> polygonPoints);
bool left(point2D a, point2D b, point2D c);
bool properIntersect(point2D a, point2D b, point2D c, point2D d);
bool between(point2D a, point2D b, point2D c);
bool improperIntersect(point2D a, point2D b, point2D c, point2D d);
bool lineProperlyIntersectsPolygon (lineSegment2D guardVertexEdge, vector<lineSegment2D> polygonEdges);
double segmentLength(point2D point1, point2D point2D);
bool lineDirection(point2D guardPoint, point2D polygonPoint);
bool polyEdgesSimple(lineSegment2D polygonEdge1, lineSegment2D polygonEdge2);
vector<lineSegment2D> extendedLines(vector<point2D> polygonPoints, point2D guardPoint);
bool isSimple(vector<point2D> polygonPoints);
bool radialSort(const point2D& pointOne, const point2D& pointTwo);
double computedAngle(point2D referencePoint, point2D newPoint);
vector<point2D> polygonLinesFunc(vector<point2D> visiblePolygonPoints);
lineSegment2D makeLine(point2D guardPoint, point2D polygonPoint);
vector<eventPoint> returnEvents(vector<point2D> polygonPoints, point2D guardPoint);
vector<point2D> visiblePolygon(vector<point2D> polygonPoints, point2D guardPoint);
vector<lineSegment2D> extendedLines(vector<point2D> polygonPoints, point2D guardPoint);
vector<vector<point2D> > triangulatedPolygon(vector<point2D> polygonPoints, point2D guardPoint);

//add any functions you might need to operate on these basic types

#endif
