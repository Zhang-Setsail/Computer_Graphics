#include "transforms.h"

#include "CGL/matrix3x3.h"
#include "CGL/vector2D.h"
#include "CGL/vector3D.h"

namespace CGL {

Vector2D operator*(const Matrix3x3 &m, const Vector2D &v) {
	Vector3D mv = m * Vector3D(v.x, v.y, 1);
	return Vector2D(mv.x / mv.z, mv.y / mv.z);
}

Matrix3x3 translate(float dx, float dy) {
	// Part 3: Fill this in.
	Matrix3x3 trans_matrix = Matrix3x3::identity();
	trans_matrix(0,2) = dx;
	trans_matrix(1,2) = dy;

	return trans_matrix;
}

Matrix3x3 scale(float sx, float sy) {
	// Part 3: Fill this in.
	Matrix3x3 scale_matrix = Matrix3x3::identity();
	scale_matrix(0,0) = sx;
	scale_matrix(1,1) = sy;

	return scale_matrix;
}

// The input argument is in degrees counterclockwise
Matrix3x3 rotate(float deg) {
	// Part 3: Fill this in.
	Matrix3x3 scale_matrix = Matrix3x3::identity();
	float cos_deg = cos(deg*PI/180.0f);
	float sin_deg = sin(deg*PI/180.0f);
	scale_matrix(0,0) = cos_deg;
	scale_matrix(0,1) = -sin_deg;
	scale_matrix(1,0) = sin_deg;
	scale_matrix(1,1) = cos_deg;

	return scale_matrix;
}

}
