#include "stubs.h"

double Test_RaySphereIntersect(glm::vec3 const& P0, glm::vec3 const& V0, glm::mat4 const& T) {
	glm::mat4 Tinv = glm::inverse(T); //use Tinv to put P0 and V0 in sphere space
	glm::vec4 P0s = Tinv * glm::vec4(P0,1); //P0s = [Tinv]*[P0]
	glm::vec4 V0s = Tinv * glm::vec4(V0,0); //V0s = [Tinv]*[V0]

	/*Set up some tests for the quadratic equation:
		t = -2*P0s*V0s +- [(2*P0s*V0s)^2 - 4(V0s*V0s)(P0s*P0s)]^.5
			---------------------------------------------------
								2(V0s*V0s)
		First determine if solutions can exist:
			(2*P0s*V0s)^2 - 4(V0s*V0s)(P0s*P0s) must be positive
			else there is no solution (no intersection)

			V0s*V0s must not be zero
			else there is no solution (and input is probably bad)
	*/
	float delta = .0000001;

	float qe1 = -2.0 * glm::dot(P0s,V0s);
	float qe2 = (2.0 * glm::dot(P0s,V0s)) * (2.0 * glm::dot(P0s,V0s));
	float qe3 = 4 * glm::dot(V0s,V0s) * glm::dot(P0s,P0s);
	float qe4 = 2.0 * glm::dot(V0s,V0s);

	if(qe2 - qe3 < 0) //if square root of a negativen number
	{
		return -1; //there must be no solution
	}

	if((qe4 > -delta) && (qe4 < delta)) //if the denominator is 0
	{
		return -1; //there must be no solution
	}

	float tp = (qe1 + glm::sqrt(qe2 - qe3))/qe4;
	float tn = (qe1 - glm::sqrt(qe2 - qe3))/qe4;

	if(tp > 0 && tn > 0)
	{
		return glm::min(tn,tp);
	}
	else if(tp < 0 && tn < 0)
	{
		return -1; //this case shouldn't happen, but if it does, sphere is behind P0
	}
	else 
	{
		return glm::max(tn,tp);
	}
}

double Test_RayPolyIntersect(glm::vec3 const& P0, glm::vec3 const& V0, glm::vec3 const& p1, glm::vec3 const& p2, glm::vec3 const& p3, glm::mat4 const& T) {
	glm::mat4 Tinv = glm::inverse(T);
	glm::vec4 P0t = Tinv * glm::vec4(P0,1);
	glm::vec4 V0t = Tinv * glm::vec4(V0,0);

	glm::vec3 p12 = p2 - p1;
	glm::vec3 p23 = p3 - p2;
	glm::vec4 N = glm::vec4(glm::cross(p12,p23),0);

	/* Through algebra magic:

		t = Nx*(P0tx - P1tx) + Ny*(P0ty - P1ty) + Nz*(P0tz - P1tz)
			-----------------------------------------------
					Nx*V0tx + Ny*V0ty + Nz*V0tz

		Intuitively, and mathematically:
			Nx*V0tx + Ny*V0ty + Nz*V0tz = 0
		Represents orthagonality to the triangle's normal (of the plane):
		We'll either never hit the plane, or we are IN the plane
	*/
	float delta = .0000001;
	float denom = glm::dot(N[0],V0t[0]) + glm::dot(N[1],V0t[1]) + glm::dot(N[2],V0t[2]);
	
	if(denom > -delta && denom < delta)
	{
		return -1; //denominator = 0, t is undefined, V0 is orthogonal
	}

	/*
		If we didn't trip that return, the we know that our ray intersect the plane of the
		triangle, but not whether the intersection is inside of the triangle. To do this,
		we're going to flatten the point of intersection and triangle to a plane, and determine
		whether the (u,v)point is inside the (u,v) triangle. The simples way to do this is drop
		the axis with the largest normal magnitude. I,e, N = [0,5,3], flatten to the XZ-plane.
	*/
	glm::vec4 diff = P0t - glm::vec4(p1,1);
	float t = (glm::dot(N[0],diff[0]) + glm::dot(N[1],diff[1]) + glm::dot(N[2],diff[2]))/denom; //calculating just in case
	glm::vec4 intersectionPoint = P0t + (t * V0t);

	if(N[0] > N[1] && N[0] > N[2])
	{
		glm::vec2 ip = glm::vec2(intersectionPoint[1],intersectionPoint[2]);
		glm::vec2 tri1 = glm::vec2(p1[1],p1[2]);
		glm::vec2 tri2 = glm::vec2(p2[1],p2[2]);
		glm::vec2 tri3 = glm::vec2(p3[1],p3[2]);
		if(isInTri(ip,tri1,tri2,tri3))
		{
			return t;
		}
		else
		{
			return -1;
		}
	}
	else if (N[1] > N[0] && N[1] > N[2])
	{
		glm::vec2 ip = glm::vec2(intersectionPoint[0],intersectionPoint[2]);
		glm::vec2 tri1 = glm::vec2(p1[0],p1[2]);
		glm::vec2 tri2 = glm::vec2(p2[0],p2[2]);
		glm::vec2 tri3 = glm::vec2(p3[0],p3[2]);
		if(isInTri(ip,tri1,tri2,tri3))
		{
			return t;
		}
		else
		{
			return -1;
		}
	}
	else //N[2] > N[1] && N[2] > N[0]     or    some equalities make the difference irrelevant
	{
		glm::vec2 ip = glm::vec2(intersectionPoint[1],intersectionPoint[0]);
		glm::vec2 tri1 = glm::vec2(p1[1],p1[0]);
		glm::vec2 tri2 = glm::vec2(p2[1],p2[0]);
		glm::vec2 tri3 = glm::vec2(p3[1],p3[0]);
		if(isInTri(ip,tri1,tri2,tri3))
		{
			return t;
		}
		else
		{
			return -1;
		}
	}


}

double Test_RayCubeIntersect(glm::vec3 const& P0, glm::vec3 const& V0, glm::mat4 const& T) {
	//Lets just split our cube up into triangles, and send transformed rays + tris to
	//the triangle intersection method.
	double t = -1;
	double d = Test_RayPolyIntersect(P0, V0, glm::vec3(-.5,-.5,.5), glm::vec3(.5,-.5,.5), glm::vec3(-.5,.5,.5), T);
	if(d > 0 && d < t) {t = d;}
	d = Test_RayPolyIntersect(P0, V0, glm::vec3(.5,-.5,.5), glm::vec3(.5,.5,.5), glm::vec3(-.5,.5,.5), T);
	if(d > 0 && d < t) {t = d;}
	d = Test_RayPolyIntersect(P0, V0, glm::vec3(-.5,-.5,-.5), glm::vec3(.5,-.5,-.5), glm::vec3(-.5,.5,-.5), T);
	if(d > 0 && d < t) {t = d;}
	d = Test_RayPolyIntersect(P0, V0, glm::vec3(.5,-.5,-.5), glm::vec3(.5,.5,-.5), glm::vec3(-.5,.5,-.5), T);
	if(d > 0 && d < t) {t = d;}
	d = Test_RayPolyIntersect(P0, V0, glm::vec3(.5,-.5,.5), glm::vec3(.5,-.5,-.5), glm::vec3(.5,.5,-.5), T);
	if(d > 0 && d < t) {t = d;}
	d = Test_RayPolyIntersect(P0, V0, glm::vec3(.5,.5,-.5), glm::vec3(.5,.5,.5), glm::vec3(.5,-.5,.5), T);
	if(d > 0 && d < t) {t = d;}
	d = Test_RayPolyIntersect(P0, V0, glm::vec3(-.5,-.5,.5), glm::vec3(-.5,-.5,-.5), glm::vec3(-.5,.5,-.5), T);
	if(d > 0 && d < t) {t = d;}
	d = Test_RayPolyIntersect(P0, V0, glm::vec3(-.5,.5,-.5), glm::vec3(-.5,.5,.5), glm::vec3(-.5,-.5,.5), T);
	if(d > 0 && d < t) {t = d;}
	d = Test_RayPolyIntersect(P0, V0, glm::vec3(-.5,.5,.5), glm::vec3(-.5,.5,-.5), glm::vec3(.5,.5,-.5), T);
	if(d > 0 && d < t) {t = d;}
	d = Test_RayPolyIntersect(P0, V0, glm::vec3(-.5,.5,.5), glm::vec3(.5,.5,-.5), glm::vec3(.5,.5,.5), T);
	if(d > 0 && d < t) {t = d;}
	d = Test_RayPolyIntersect(P0, V0, glm::vec3(-.5,-.5,.5), glm::vec3(-.5,-.5,-.5), glm::vec3(.5,-.5,-.5), T);
	if(d > 0 && d < t) {t = d;}
	d = Test_RayPolyIntersect(P0, V0, glm::vec3(-.5,-.5,.5), glm::vec3(.5,-.5,-.5), glm::vec3(.5,-.5,.5), T);
	if(d > 0 && d < t) {t = d;}

	//Hurts to read. I'll use a more efficient algorithm later.

	return t;
}

double Test_RayCylinderIntersect(glm::vec3 const& P0, glm::vec3 const& V0, glm::mat4 const& T) {
	// TODO fill this in so it calls your own ray-casting function.
	// See the documentation of this function in stubs.h.

	return -1;
}

bool isInTri(glm::vec2 p,glm::vec2 p1,glm::vec2 p2,glm::vec2 p3)
{
	glm::vec3 p12 = glm::vec3(p2 - p1, 0);
	glm::vec3 p13 = glm::vec3(p3 - p1, 0);
	glm::vec3 p1p = glm::vec3(p - p1, 0);
	glm::vec3 p1cross23 = glm::cross(p12,p13);
	glm::vec3 p1cross2p = glm::cross(p12,p1p);
	bool cross1 = (p1cross23[2] < 0 && p1cross2p[2] < 0) || (p1cross23[2] > 0 && p1cross2p[2] > 0);

	glm::vec3 p23 = glm::vec3(p3 - p2, 0);
	glm::vec3 p21 = glm::vec3(p1 - p2, 0);
	glm::vec3 p2p = glm::vec3(p - p2, 0);
	glm::vec3 p2cross31 = glm::cross(p23,p21);
	glm::vec3 p2cross3p = glm::cross(p23,p2p);
	bool cross2 = (p2cross31[2] < 0 && p2cross3p[2] < 0) || (p2cross31[2] > 0 && p2cross3p[2] > 0);

	glm::vec3 p31 = glm::vec3(p3 - p1, 0);
	glm::vec3 p32 = glm::vec3(p2 - p1, 0);
	glm::vec3 p3p = glm::vec3(p - p3, 0);
	glm::vec3 p3cross12 = glm::cross(p31,p32);
	glm::vec3 p3cross1p = glm::cross(p31,p3p);
	bool cross3 = (p3cross12[2] < 0 && p3cross1p[2] < 0) || (p3cross12[2] > 0 && p3cross1p[2] > 0);

	return (cross1 && cross2 && cross3);
}