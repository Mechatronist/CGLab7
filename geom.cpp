#include "geom.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>

#include <vector>
#include <map>
#include <utility>
#include <math.h>

using namespace std; 

int WINDOWSIZE = 750;

bool isSimple(vector<point2D> polygonPoints) {
  
	vector<lineSegment2D> polygonEdges = makePolygonEdges(polygonPoints);

	for (int i = 0; i < polygonEdges.size(); i++) {
		for (int j = 0; j < polygonEdges.size(); j++) {

			if (properIntersect(polygonEdges[i].p1, polygonEdges[i].p2, polygonEdges[j].p1, polygonEdges[j].p2)) {
		 		return false;
		 	}
		}
	} 

	return true;
}

vector<point2D> polygonShift(vector<point2D> originalPolygonPoints, point2D guardPoint) {

	for (int i = 0; i < originalPolygonPoints.size(); i++) {
		originalPolygonPoints[i].x -= guardPoint.x;
		originalPolygonPoints[i].y -= guardPoint.y;
	} 

	return originalPolygonPoints;

}

double computeXIntersection (int firstPointX, int firstPointY, int secondPointX, int secondPointY) {
	return ((firstPointX * secondPointY - secondPointX * firstPointY)) / ((double)(secondPointY - firstPointY));

}

bool isInPolygon(vector<point2D> polygonPoints, point2D guardPoint) { //Need last mouse click for point

	int previousPoint; // = i1
	int numberOfPoints = polygonPoints.size();
	double xAxisIntersection;
	int Rcross = 0;
	int Lcross = 0;
	bool Rstrad, Lstrad;

	vector<point2D> shiftedPoints;

	shiftedPoints = polygonShift(polygonPoints, guardPoint); //Points are shifted by guardPoint

	for (int i = 0; i < polygonPoints.size(); i++) { // i = first edge point

		if (shiftedPoints[i].x == 0 && shiftedPoints[i].y == 0) {
			return false;
		}

		previousPoint = (i + numberOfPoints - 1) % numberOfPoints;

		Rstrad = (shiftedPoints[i].y > 0) != (shiftedPoints[previousPoint].y > 0);
		Lstrad = (shiftedPoints[i].y < 0) != (shiftedPoints[previousPoint].y < 0);

		if (Rstrad || Lstrad) {

			xAxisIntersection = computeXIntersection(shiftedPoints[i].x, shiftedPoints[i].y, shiftedPoints[previousPoint].x, shiftedPoints[previousPoint].y);

			if (Rstrad && xAxisIntersection > 0) {
				Rcross += 1;
			}

			if (Lstrad && xAxisIntersection > 0) {
				Lcross += 1;
			}

		}

	}

	if ((Rcross % 2) != (Lcross % 2)) {
		return false;
	}

	if ( (Rcross % 2) == 1) {
		return true;
	}

	else {
		return false;
	}

}

bool collinear(point2D a, point2D b, point2D c) {
	return (a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y)) == 0;
}

int Area2(point2D a, point2D b, point2D c) {
	return (b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y);
}

bool Xor(bool x, bool y) {
	return !x ^ !y;
}

bool left(point2D a, point2D b, point2D c) {
	return Area2(a,b,c) > 0;
}

bool properIntersect(point2D a, point2D b, point2D c, point2D d) {
	if (collinear(a,b,c) || collinear(a,b,d) || collinear(c,d,a) || collinear(c,d,b)) {
		return false;
	}

	return Xor(left(a,b,c), left(a,b,d)) && Xor(left(c,d,a), left(c,d,b));
}

bool between(point2D a, point2D b, point2D c) {
	
	if (!collinear(a,b,c)) {
		return false;
	}

	if (a.x != b.x) {
		return ((a.x <= c.x) && (c.x <= b.x)) || ((a.x >= c.x) && (c.x >= b.x));
	}

	else {
		return ((a.y <= c.y) && (c.y <= b.y)) || ((a.y >= c.y) && (c.y >= b.y));
	}

}

bool improperIntersect(point2D a, point2D b, point2D c, point2D d) {
	return between(a, b, c) || between(a, b, d) || between(c, d, a) || between(c, d, b);
}

double segmentLength(point2D point1, point2D point2) {
	return sqrt((pow((point1.x - point2.x), 2)) + (pow((point1.y - point2.y), 2)));
}

lineSegment2D extendLine(point2D polygonPoint, point2D guardPoint) {

	lineSegment2D extendedLine;
	point2D newEndPoint;

	newEndPoint.x = WINDOWSIZE*(polygonPoint.x - guardPoint.x) + polygonPoint.x;
	newEndPoint.y = WINDOWSIZE*(polygonPoint.y - guardPoint.y) + polygonPoint.y;

	extendedLine.p1 = guardPoint;
	extendedLine.p2 = newEndPoint;

	return extendedLine;

}

lineSegment2D extendLineBothDirections(point2D polygonPoint, point2D guardPoint) {

	lineSegment2D extendedLine;
	point2D newEndPoint, newEndPoint2;

	newEndPoint.x = WINDOWSIZE*(polygonPoint.x - guardPoint.x) + polygonPoint.x;
	newEndPoint.y = WINDOWSIZE*(polygonPoint.y - guardPoint.y) + polygonPoint.y;

	newEndPoint2.x = WINDOWSIZE*(guardPoint.x - polygonPoint.x) + guardPoint.x;
	newEndPoint2.y = WINDOWSIZE*(guardPoint.y - polygonPoint.y) + guardPoint.y;

	extendedLine.p1 = newEndPoint;
	extendedLine.p2 = newEndPoint2;

	return extendedLine;

}

bool linesIntersect(lineSegment2D guardVertexEdge, lineSegment2D polygonEdge) {

	if (properIntersect(guardVertexEdge.p1, guardVertexEdge.p2, polygonEdge.p1, polygonEdge.p2)) {
		return true;
	}

	else if (improperIntersect(guardVertexEdge.p1, guardVertexEdge.p2, polygonEdge.p1, polygonEdge.p2)) {
		return true;
	}

	else {
		return false;
	}

}

double slopeOfLine(lineSegment2D line) {
	return (line.p1.y - line.p2.y) / (line.p1.x - line.p2.x);
}

double yIntercept(lineSegment2D line, double slope) {
	return line.p1.y - (slope * line.p1.x);
}

bool grazePoint(point2D polygonPointV1, point2D polygonPointV2, lineSegment2D guardVertexEdge) {

	lineSegment2D prevNextLine;

	prevNextLine.p1 = polygonPointV1;
	prevNextLine.p2 = polygonPointV2;

	if (properIntersect(prevNextLine.p1, prevNextLine.p2, guardVertexEdge.p1, guardVertexEdge.p2)) {
		return false;
	}

	else {
		return true;
	}

}

point2D intersectionPoint(point2D a, point2D b, point2D c, point2D d) {

	double s;
	double num, denom;
	point2D intersectionPoint;

	denom = a.x * (double) (d.y - c.y) + b.x * (double) (c.y - d.y) + d.x * (double) (b.y - a.y) + c.x * (double) (a.y - b.y);

	num = a.x * (double) (d.y - c.y) + c.x * (double) (a.y - d.y) + d.x * (double) (c.y - a.y);

	s = num / denom;

	intersectionPoint.x = a.x + s * (b.x - a.x);
	intersectionPoint.y = a.y + s * (b.y - a.y);

	return intersectionPoint;

}

bool eventSorterDistance(const eventPoint& pointOne, const eventPoint& pointTwo) {
	return pointOne.distanceFromGuard < pointTwo.distanceFromGuard;
}

bool pointSorterDistance(const point2D& pointOne, const point2D& pointTwo) {
	return pointOne.distanceFromGuard < pointTwo.distanceFromGuard;
}

bool pointSorterNegativeDistance(const point2D& pointOne, const point2D& pointTwo) {
	return pointOne.distanceFromGuard > pointTwo.distanceFromGuard;
}

double pointDistance(point2D event, point2D guardPoint) {
	return sqrt((pow((event.x - guardPoint.x), 2)) + (pow((event.y - guardPoint.y), 2)));
}

vector<eventPoint> returnEvents(vector<point2D> polygonPoints, point2D guardPoint) { //Returns polygon visible to guard.

	vector<point2D> visiblePolygon;
	vector<lineSegment2D> polygonEdges = makePolygonEdges(polygonPoints);
	vector<eventPoint> intersectionEvents;

	lineSegment2D guardVertexEdge;
	eventPoint newEvent;

	for (int i = 0; i < polygonPoints.size(); i++) { //Need to check if point is on a vertex. If so, just add it.

		guardVertexEdge = extendLine(polygonPoints[i], guardPoint);

		for (int j = 0; j < polygonEdges.size(); j++) {

			if (linesIntersect(guardVertexEdge, polygonEdges[j]) == true) {
		 		newEvent.point = intersectionPoint(guardVertexEdge.p1, guardVertexEdge.p2, polygonEdges[j].p1, polygonEdges[j].p2);
		 		newEvent.intersectionType = properIntersect(guardVertexEdge.p1, guardVertexEdge.p2, polygonEdges[j].p1, polygonEdges[j].p2);
		 		newEvent.distanceFromGuard = pointDistance(newEvent.point, guardPoint);
		 		intersectionEvents.push_back(newEvent);

		 	}

		}

	}
	return intersectionEvents;

}

bool pointExists(vector<point2D> visiblePolygon, point2D newPoint) {

	for (int i = 0; i < visiblePolygon.size(); i++) {
		if (visiblePolygon[i].x == newPoint.x && visiblePolygon[i].y == newPoint.y) {
			return true;
		}
	}

	return false;

}

bool eventExists(vector<eventPoint> events, point2D newPoint) {

	for (int i = 0; i < events.size(); i++) {
		if (events[i].point.x == newPoint.x && events[i].point.y == newPoint.y) {
			return true;
		}
	}

	return false;

}

vector<lineSegment2D> makePolygonEdges(vector<point2D> polygonPoints) { //Returns all edges in polygon

	vector<lineSegment2D> polygonEdges;
	lineSegment2D newEdge;

	for (int i = 0; i < polygonPoints.size() - 1; i++) {
		newEdge.p1 = polygonPoints[i];
		newEdge.p2 = polygonPoints[i+1];
		polygonEdges.push_back(newEdge);
	}

	newEdge.p1 = polygonPoints.back();
	newEdge.p2 = polygonPoints.front();
	polygonEdges.push_back(newEdge);

	return polygonEdges;

}

double computedAngle(point2D referencePoint, point2D newPoint) {

	double computeAngle = atan2((newPoint.y - referencePoint.y), (newPoint.x - referencePoint.x));

	if (computeAngle < 0) {
        computeAngle = computeAngle + 2*M_PI;
    }

    return computeAngle;

}

bool pointComparer(const point2D& pointOne, const point2D& pointTwo) {

  if (pointOne.x == pointTwo.x && pointOne.y == pointTwo.y) {
    return true;
  }

  return false;

}

bool distanceSort(const point2D& pointOne, const point2D& pointTwo) {
  return pointOne.distanceFromPreviousPoint < pointTwo.distanceFromPreviousPoint;

}

void remove_coincident_points(vector<point2D> *points) {
  vector<point2D>::iterator unique_end = (unique((*points).begin(), (*points).end(), pointComparer));
  (*points).erase(unique_end, (*points).end());
  
}

vector<point2D> reformPolygon(vector<point2D> visiblePolygonPoints, vector<lineSegment2D> polygonEdges) {

	vector<point2D> edgePoints, reformedPolygon;

	for (int i = 0; i < polygonEdges.size(); i++) {  //Ordered polygon edges.
		for (int j = 0; j < visiblePolygonPoints.size(); j++) {

			if (visiblePolygonPoints[j].intersectionEdge == i) {
				edgePoints.push_back(visiblePolygonPoints[j]);
			}

		}

		if (edgePoints.size() == 0) { //No visible points on that edge.
			continue;
		}

		else if (edgePoints.size() == 1) { //Only one visible point on that edge.
			reformedPolygon.push_back(edgePoints[0]);
		}

		else { 

			sort(edgePoints.begin(), edgePoints.end(), distanceSort); 

			for (int j = 0; j < edgePoints.size(); j++) {
				reformedPolygon.push_back(edgePoints[j]);

			}

		}

		edgePoints.clear();

	}

	return reformedPolygon;

}

void initializeNewEvent(eventPoint *newEvent, int intersectionEdge, lineSegment2D polygonEdge, lineSegment2D guardVertexEdge, point2D guardPoint) {
	(*newEvent).intersectionEdge = intersectionEdge;
	(*newEvent).distanceFromPreviousPoint = pointDistance((*newEvent).point, polygonEdge.p1);
	(*newEvent).intersectionType = properIntersect(guardVertexEdge.p1, guardVertexEdge.p2, polygonEdge.p1, polygonEdge.p2);
	(*newEvent).distanceFromGuard = pointDistance((*newEvent).point, guardPoint);

}

void initializePoint(eventPoint *newEvent, point2D guardPoint, point2D prevPolygonPoint) {
	(*newEvent).point.distanceFromGuard = (*newEvent).distanceFromGuard;
	(*newEvent).point.guardAngle = computedAngle(guardPoint, (*newEvent).point);
	(*newEvent).point.intersectionType = (*newEvent).intersectionType;
	(*newEvent).point.left = left(guardPoint, (*newEvent).point, prevPolygonPoint);
	(*newEvent).point.intersectionEdge = (*newEvent).intersectionEdge;
	(*newEvent).point.distanceFromPreviousPoint = (*newEvent).distanceFromPreviousPoint;
}

vector<point2D> visiblePolygon(vector<point2D> polygonPoints, point2D guardPoint) { //Returns polygon visible to guard.

	vector<point2D> visiblePolygon;
	vector<lineSegment2D> polygonEdges = makePolygonEdges(polygonPoints);
	vector<eventPoint> intersectionEvents, tempContainer;

	lineSegment2D guardVertexEdge;
	eventPoint newEvent;

	int prevPointIndex, nextPointIndex;

	remove_coincident_points(&polygonPoints);

	for (int i = 0; i < polygonPoints.size(); i++) { //Need to check if point is on a vertex. If so, just add it.

		if (i == 0) {
			prevPointIndex = polygonPoints.size() - 1;
			nextPointIndex = 1;
		}

		else if (i == polygonPoints.size() - 1) {
			prevPointIndex = polygonPoints.size() - 2;
			nextPointIndex = 0;
		}

		else {
			prevPointIndex = i - 1;
			nextPointIndex = i + 1;
		}

		guardVertexEdge = extendLine(polygonPoints[i], guardPoint);

		for (int j = 0; j < polygonEdges.size(); j++) {

			if (linesIntersect(guardVertexEdge, polygonEdges[j])) {
		 		newEvent.point = intersectionPoint(guardVertexEdge.p1, guardVertexEdge.p2, polygonEdges[j].p1, polygonEdges[j].p2);

		 		if (eventExists(intersectionEvents, newEvent.point)) {
					continue;
				}

				initializeNewEvent(&newEvent, j, polygonEdges[j], guardVertexEdge, guardPoint);
		 		intersectionEvents.push_back(newEvent);

		 	}

		}

		sort(intersectionEvents.begin(), intersectionEvents.end(), eventSorterDistance);

		for (int k = 0; k < intersectionEvents.size(); k++) {

			initializePoint(&intersectionEvents[k], guardPoint, polygonPoints[prevPointIndex]);
			guardVertexEdge = extendLineBothDirections(polygonPoints[i], guardPoint);

			if (intersectionEvents[k].intersectionType == 0 && grazePoint(polygonPoints[prevPointIndex], polygonPoints[nextPointIndex], guardVertexEdge)) {
				if (!pointExists(visiblePolygon, intersectionEvents[k].point)) {
					visiblePolygon.push_back(intersectionEvents[k].point);
				}
			}

			else {  //IntersectionType is proper or improper and end
				if (!pointExists(visiblePolygon, intersectionEvents[k].point)) {
					visiblePolygon.push_back(intersectionEvents[k].point);
				}

				break;
			}

		}

		intersectionEvents.clear();

	} //Returns polygon visible to guard.

	visiblePolygon = reformPolygon(visiblePolygon, polygonEdges);

	return visiblePolygon;
}


vector<vector<point2D> > triangulatedPolygon(vector<point2D> polygonPoints, point2D guardPoint) {

	vector<vector<point2D> > triangulatedPolygon;
	vector<point2D> visiblePolygonPoints, singleTriangle;
	point2D p1, p2;

	visiblePolygonPoints = visiblePolygon(polygonPoints, guardPoint);

	for (int i = 0; i < visiblePolygonPoints.size() - 1; i ++) {
		p1 = visiblePolygonPoints[i];
		p2 = visiblePolygonPoints[i+1];

		singleTriangle.push_back(p1);
		singleTriangle.push_back(p2);
		singleTriangle.push_back(guardPoint);

		triangulatedPolygon.push_back(singleTriangle);

		singleTriangle.clear();

	}

	p1 = visiblePolygonPoints.back();
	p2 = visiblePolygonPoints.front();

	singleTriangle.push_back(p1);
	singleTriangle.push_back(guardPoint);
	singleTriangle.push_back(p2);

	triangulatedPolygon.push_back(singleTriangle);

	return triangulatedPolygon;

}