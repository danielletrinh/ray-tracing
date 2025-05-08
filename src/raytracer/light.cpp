#include "light.h"

#include <sstream>
#include <string>
#include <cmath>
#include <cstdlib>

#define ran() ( double( rand() % 32768 ) / 32768 )

Light::Light() {
    sample = rand();
    next = nullptr;
    lightPrimitive = nullptr;
}

void Light::Input( std::string var , std::stringstream& fin ) {
    if ( var == "color=" ) color.Input( fin );
}

// ==========================================

void PointLight::Input( std::string var , std::stringstream& fin ) {
    if ( var == "O=" ) O.Input( fin );
    Light::Input( var , fin );
}

double PointLight::CalnShade( Vector3 C , Primitive* primitive_head , int shade_quality ) {
    /* For point light, shade_quality is of no use: we don't need to sample. */

    // light ray from diffuse point to light source
    Vector3 V = O - C;
    double dist = V.Module();

    // if light ray collide any object, light source produce no shade to diffuse light
    for (Primitive* now = primitive_head ; now != NULL ; now = now->GetNext()) {
        CollidePrimitive tmp = now->Collide(C, V);
        if (dist - tmp.dist > EPS) {
            return 0;
        }
    }

    return 1;
}

// ==========================================

void SquareLight::Input( std::string var , std::stringstream& fin ) {
    if ( var == "O=" ) O.Input( fin );
    if ( var == "Dx=" ) Dx.Input( fin );
    if ( var == "Dy=" ) Dy.Input( fin );
    Light::Input( var , fin );
}

double SquareLight::CalnShade( Vector3 C , Primitive* primitive_head , int shade_quality ) {
    int shade = 0;

    for (int i = 0; i < shade_quality; i++) {
        // sample a point light from light primitive
        Vector3 randO = GetRandPointLight(C);

        // light ray from diffuse point to point light
        Vector3 V = randO - C;
        double dist = V.Module();

        int addShade = 1;
        // if light ray collide any object before reaching the light, point light produce no shade to diffuse light
        for (Primitive* now = primitive_head; now != NULL; now = now->GetNext()) {
            // don't collide with myself!
            // this shouldn't have been a trouble, but there is COMPUTATIONAL ERROR!
            if (now == lightPrimitive) {
                continue;
            }
            CollidePrimitive tmp = now->Collide(C, V);
            if (dist - tmp.dist > EPS) {
                addShade = 0;
                break;
            }
        }
        shade += addShade;
    }
    return double(shade) / shade_quality;
}

Primitive* SquareLight::CreateLightPrimitive()
{
    PlaneAreaLightPrimitive* res = new PlaneAreaLightPrimitive(O, Dx, Dy, color);
    lightPrimitive = res;
    return res;
}

Vector3 SquareLight::GetRandPointLight(const Vector3& crashPoint)
{
    return O + Dx * (2 * ran() - 1) + Dy * (2 * ran() - 1);
}


// ==========================================

void SphereLight::Input( std::string var , std::stringstream& fin ) {
     if ( var == "O=" ) O.Input( fin );
     if ( var == "R=" ) fin>>R;
     Light::Input( var , fin );
}

double SphereLight::CalnShade( Vector3 C , Primitive* primitive_head , int shade_quality ) {
    int shade = 0;

    for (int i = 0; i < shade_quality; i++) {
        Vector3 randO = GetRandPointLight(C);
        Vector3 V = randO - C;
        double dist = V.Module();

        bool visible = true;
        for (Primitive* now = primitive_head; now != nullptr; now = now->GetNext()) {
            if (now == lightPrimitive) continue;
            CollidePrimitive tmp = now->Collide(C, V);
            if (tmp.dist < dist - EPS) {
                visible = false;
                break;
            }
        }
        shade += visible ? 1 : 0;
    }
    return double(shade) / shade_quality;
}

Primitive* SphereLight::CreateLightPrimitive()
{
    SphereLightPrimitive* res = new SphereLightPrimitive(O, R, color);
    lightPrimitive = res;
    return res;
}

Vector3 SphereLight::GetRandPointLight(const Vector3& crashPoint)
{
    double theta = 2 * M_PI * ran();
    double z = 2 * ran() - 1;  // Uniform in [-1, 1]
    double r = sqrt(1 - z*z);  // Radius in xy-plane

    Vector3 dir = Vector3(r * cos(theta), r * sin(theta), z);
    return O + R * dir;
}