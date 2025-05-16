#ifndef PRIMITIVE_H
#define PRIMITIVE_H

#include<iostream>
#include<sstream>
#include<string>
#include<vector>

#include"../helpers/color.h"
#include"../helpers/vector3.h"
#include"../helpers/bmp.h"

extern const double EPS;
extern const double PI;
const double BIG_DIST = 1e10;

class Blur {
public:
    virtual std::pair<double, double> GetXY() = 0;
};

class ExpBlur : public Blur {
public:
    /* return x & y coord of a random point inside a unit sphere,
       with radius following an exponential sampling.
    */
    std::pair<double, double> GetXY();
};

class Material {
public:
    Color color , absor;
    double refl , refr;
    double diff , spec;
    double rindex;
    double drefl;
    Bmp* texture;
    Blur* blur;

    Material();
    ~Material() {}

    void Input( std::string , std::stringstream& );
};

struct CollidePrimitive;

class Primitive {
protected:
    int sample;
    Material* material;
    Primitive* next;

public:

    Primitive();
    Primitive( const Primitive& );
    ~Primitive();

    int GetSample() { return sample; }
    Material* GetMaterial() { return material; }
    Primitive* GetNext() { return next; }
    void SetNext( Primitive* primitive ) { next = primitive; }

    virtual void Input( std::string , std::stringstream& );
    virtual CollidePrimitive Collide( Vector3 ray_O , Vector3 ray_V ) = 0;
    virtual Color GetTexture(Vector3 crash_C) = 0;
    virtual bool IsLightPrimitive() const { return false; }
};

struct CollidePrimitive {
    bool isCollide;
    Primitive* collide_primitive;
    Vector3 N , C;
    double dist;
    bool front;
    CollidePrimitive(){isCollide = false; collide_primitive = NULL; dist = BIG_DIST;}
    Color GetTexture(){return collide_primitive->GetTexture(C);}
};

class Sphere : public Primitive {
protected:
    Vector3 O , De , Dc;
    double R;

public:
    Sphere();
    ~Sphere() {}

    void Input( std::string , std::stringstream& );
    CollidePrimitive Collide( Vector3 ray_O , Vector3 ray_V );
    Color GetTexture(Vector3 crash_C);
};

class SphereLightPrimitive : public Sphere{
public:
    SphereLightPrimitive(Vector3 pO, double pR, Color color) : Sphere()
    {O = pO; R = pR; material->color = color; }
    virtual bool IsLightPrimitive() const { return true; }
};

class Plane : public Primitive {
protected:
    Vector3 N , Dx , Dy;
    double R;

public:
    Plane() : Primitive() {}
    ~Plane() {}

    void Input( std::string , std::stringstream& );
    CollidePrimitive Collide( Vector3 ray_O , Vector3 ray_V );
    Color GetTexture(Vector3 crash_C);
};

class Square : public Primitive {
protected:
    Vector3 O , Dx , Dy;

public:
    Square() : Primitive() {}
    ~Square() {}

    void Input( std::string , std::stringstream& );
    CollidePrimitive Collide( Vector3 ray_O , Vector3 ray_V );
    Color GetTexture(Vector3 crash_C);
};

class PlaneAreaLightPrimitive : public Square {
public:
    explicit PlaneAreaLightPrimitive(Vector3 pO, Vector3 pDx, Vector3 pDy, Color color): Square()
    {O = pO; Dx = pDx; Dy = pDy; material->color = color; }
    virtual bool IsLightPrimitive() const { return true; }
};

class Cylinder : public Primitive {
    Vector3 O1, O2;
    double R;

public:
    Cylinder() : Primitive() {}
    explicit Cylinder(Vector3 pO1, Vector3 pO2, double pR) : Primitive() {O1 = pO1; O2 = pO2; R = pR; }
    ~Cylinder() {}

    void Input( std::string , std::stringstream& );
    CollidePrimitive Collide( Vector3 ray_O , Vector3 ray_V );
    Color GetTexture(Vector3 crash_C);
};

class Bezier : public Primitive {
    Vector3 O1, O2;
    Vector3 N, Nx, Ny;
    std::vector<double> R;
    std::vector<double> Z;
    int degree;
    Cylinder* boundingCylinder;

public:
    Bezier() : Primitive() {boundingCylinder = NULL; degree = -1; N = Vector3(0, 0, 1); Nx = Vector3(1, 0, 0); Ny = Vector3(0, 1, 0);}
    ~Bezier() {delete boundingCylinder;}

    void Input( std::string , std::stringstream& );
    CollidePrimitive Collide( Vector3 ray_O , Vector3 ray_V );
    Color GetTexture(Vector3 crash_C);

private:
    std::pair<double, double> valueAt(double u);
    std::pair<double, double> valueAt(double u, const std::vector<double>& xs, const std::vector<double>& ys);
};

class Voxel : public Primitive {
protected:
    Vector3 O;          // Origin (bottom-left-front corner)
    Vector3 size;       // Size of each voxel
    int nx, ny, nz;     // Grid dimensions
    std::vector<bool> grid; // 3D grid stored as 1D array

public:
    Voxel() : Primitive(), nx(0), ny(0), nz(0) {}
    ~Voxel() {grid.clear();}

    void Input(std::string var, std::stringstream& fin);
    CollidePrimitive Collide(Vector3 ray_O, Vector3 ray_V);
    Color GetTexture(Vector3 crash_C);

private:
    int GetIndex(int x, int y, int z) const {
        return x + nx * (y + ny * z);
    }
    bool IsVoxelSet(int x, int y, int z) const {
        if (x < 0 || x >= nx || y < 0 || y >= ny || z < 0 || z >= nz) return false;
        return grid[GetIndex(x, y, z)];
    }
};

#endif