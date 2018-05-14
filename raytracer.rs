// #if defined __linux__ || defined __APPLE__
// // "Compiled for Linux
// #else
// // Windows doesn't define these values by default, Linux does
// #define M_PI 3.141592653589793
// #define INFINITY 1e8
// #endif
// #define MAX_RAY_DEPTH 5

use std::ops::Add;
use std::ops::Mul;
use std::ops::Neg;
use std::ops::Sub;

const PI: f32 = 3.141592653589793;
const INFINITY: f32 = 1e8;
const MAX_RAY_DEPTH: u64 = 5;

#[derive(Clone, Copy)]
struct Vec3 {
    x: f32,
    y: f32,
    z: f32,
}
impl Add for Vec3 {
    type Output = Vec3;

    fn add(self, other: Vec3) -> Vec3 {
        Vec3 {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}
impl Mul for Vec3 {
    type Output = Vec3;

    fn mul(self, other: Vec3) -> Vec3 {
        Vec3 {
            x: self.x * other.x,
            y: self.y * other.y,
            z: self.z * other.z,
        }
    }
}
impl Sub for Vec3 {
    type Output = Vec3;

    fn sub(self, other: Vec3) -> Vec3 {
        Vec3 {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}
impl Neg for Vec3 {
    type Output = Vec3;

    fn neg(self) -> Vec3 {
        Vec3 {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}

// template<typename T>
// class Vec3
// {
// public:
//     T x, y, z;

impl Vec3 {
    //     Vec3() : x(T(0)), y(T(0)), z(T(0)) {}
    fn build0() -> Vec3 {
        Vec3 {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        }
    }

    //     Vec3(T xx) : x(xx), y(xx), z(xx) {}
    fn build1(v: f32) -> Vec3 {
        Vec3 { x: v, y: v, z: v }
    }
    //     Vec3(T xx, T yy, T zz) : x(xx), y(yy), z(zz) {}

    //     T length() const { return sqrt(length2()); }
    fn length(&self) -> f32 {
        self.length2().sqrt()
    }
    // length2() const { return x * x + y * y + z * z; }
    fn length2(&self) -> f32 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }
    //     Vec3& normalize()
    //     {
    //         T nor2 = length2();
    //         if (nor2 > 0) {
    //             T invNor = 1 / sqrt(nor2);
    //             x *= invNor, y *= invNor, z *= invNor;
    //         }
    //         return *this;
    //     }
    fn normalize(&mut self) {
        let nor2 = self.length2();
        if nor2 > 0.0 {
            let inv_nor = 1.0 / nor2.sqrt();
            self.x *= inv_nor;
            self.y *= inv_nor;
            self.z *= inv_nor;
        }
    }
    // T dot(const Vec3<T> &v) const { return x * v.x + y * v.y + z * v.z; }
    fn dot(&self, v: &Vec3) -> f32 {
        self.x * v.x + self.y * v.y + self.z * v.z
    }
    //     Vec3<T> operator - (const Vec3<T> &v) const { return Vec3<T>(x - v.x, y - v.y, z - v.z); }

    fn sub(&self, v: Vec3) -> Vec3 {
        Vec3 {
            x: self.x - v.x,
            y: self.y - v.y,
            z: self.z - v.z,
        }
    }

    //     Vec3<T> operator * (const T &f) const { return Vec3<T>(x * f, y * f, z * f); }
    //     Vec3<T> operator * (const Vec3<T> &v) const { return Vec3<T>(x * v.x, y * v.y, z * v.z); }
}

//     Vec3<T> operator + (const Vec3<T> &v) const { return Vec3<T>(x + v.x, y + v.y, z + v.z); }
//     Vec3<T>& operator += (const Vec3<T> &v) { x += v.x, y += v.y, z += v.z; return *this; }
//     Vec3<T>& operator *= (const Vec3<T> &v) { x *= v.x, y *= v.y, z *= v.z; return *this; }
//     Vec3<T> operator - () const { return Vec3<T>(-x, -y, -z); }
//     friend std::ostream & operator << (std::ostream &os, const Vec3<T> &v)
//     {
//         os << "[" << v.x << " " << v.y << " " << v.z << "]";
//         return os;
//     }
// };

// class Sphere
// {
// public:
//     Vec3f center;                           /// position of the sphere
//     float radius, radius2;                  /// sphere radius and radius^2
//     Vec3f surfaceColor, emissionColor;      /// surface color and emission (light)
//     float transparency, reflection;         /// surface transparency and reflectivity

#[derive(Clone, Copy)]
struct Sphere {
    center: Vec3,
    radius: f32,
    radius2: f32,
    surface_color: Vec3,
    emission_color: Vec3,
    transparency: f32,
    reflection: f32,
}

//     Sphere(
//         const Vec3f &c,
//         const float &r,
//         const Vec3f &sc,
//         const float &refl = 0,
//         const float &transp = 0,
//         const Vec3f &ec = 0) :
//         center(c), radius(r), radius2(r * r), surfaceColor(sc), emissionColor(ec),
//         transparency(transp), reflection(refl)
//     { /* empty */ }

impl Sphere {
    fn build(c: Vec3, r: f32, sc: Vec3, ec: Vec3, refl: f32, transp: f32) -> Sphere {
        Sphere {
            center: c,
            radius: r,
            radius2: r * r,
            surface_color: sc,
            emission_color: ec,
            transparency: transp,
            reflection: refl,
        }
    }

    //     //[comment]
    //     // Compute a ray-sphere intersection using the geometric solution
    //     //[/comment]
    //     bool intersect(const Vec3f &rayorig, const Vec3f &raydir, float &t0, float &t1) const
    //     {
    //         Vec3f l = center - rayorig;
    //         float tca = l.dot(raydir);
    //         if (tca < 0) return false;
    //         float d2 = l.dot(l) - tca * tca;
    //         if (d2 > radius2) return false;
    //         float thc = sqrt(radius2 - d2);
    //         t0 = tca - thc;
    //         t1 = tca + thc;

    //         return true;
    //     }
    fn intersect(&self, rayorig: Vec3, raydir: Vec3, t0: &mut f32, t1: &mut f32) -> bool {
        let l: Vec3 = self.center.sub(rayorig);
        let tca: f32 = l.dot(&raydir);
        if tca < 0.0 {
            return false;
        }
        let d2: f32 = l.dot(&l) - tca * tca;
        let thc: f32 = (self.radius - d2).sqrt();
        let t0: f32 = tca - thc;
        let t1: f32 = tca + thc;
        return true;
    }
}

// float mix(const float &a, const float &b, const float &mix)
// {
//     return b * mix + a * (1 - mix);
// }

fn mix(a: f32, b: f32, mix: f32) -> f32 {
    b * mix + a * (1.0 - mix)
}

// Vec3f trace(
//     const Vec3f &rayorig,
//     const Vec3f &raydir,
//     const std::vector<Sphere> &spheres,
//     const int &depth)

fn trace(rayorig: Vec3, raydir: Vec3, spheres: &Vec<Sphere>, depth: u64) -> Vec3 {
    // {
    //     //if (raydir.length() != 1) std::cerr << "Error " << raydir << std::endl;
    if raydir.length() != 1.0 {
        panic!("raydir length wasn't 1");
    }
    //     float tnear = INFINITY;
    //     const Sphere* sphere = NULL;
    //     // find intersection of this ray with the sphere in the scene
    //     for (unsigned i = 0; i < spheres.size(); ++i) {
    //         float t0 = INFINITY, t1 = INFINITY;
    //         if (spheres[i].intersect(rayorig, raydir, t0, t1)) {
    //             if (t0 < 0) t0 = t1;
    //             if (t0 < tnear) {
    //                 tnear = t0;
    //                 sphere = &spheres[i];
    //             }
    //         }
    //     }

    let mut tnear = INFINITY;
    let mut osphere: Option<Sphere> = None;
    for vsphere in spheres.iter() {
        let mut t0 = INFINITY;
        let mut t1 = INFINITY;
        if vsphere.intersect(rayorig, raydir, &mut t0, &mut t1) {
            if t0 < 0.0 {
                t0 = t1;
            }
            if t0 < tnear {
                tnear = t0;
                let cloned_sphere = vsphere.clone();
                let o_cs = Some(cloned_sphere);
                osphere = o_cs;
            }
        }
    }
    //     // if there's no intersection return black or background color
    //     if (!sphere) return Vec3f(2);
    let sphere = match osphere {
        None => return Vec3::build1(2.0),
        Some(sphere) => sphere,
    };
    //     Vec3f surfaceColor = 0; // color of the ray/surfaceof the object intersected by the ray
    //     Vec3f phit = rayorig + raydir * tnear; // point of intersection
    //     Vec3f nhit = phit - sphere->center; // normal at the intersection point
    //     nhit.normalize(); // normalize normal direction
    let surface_color = Vec3::build1(0.0);
    let phit = rayorig + raydir * Vec3::build1(tnear);
    let mut nhit = phit - sphere.center;
    nhit.normalize();
    //     float bias = 1e-4; // add some bias to the point from which we will be tracing
    //     bool inside = false;
    //     if (raydir.dot(nhit) > 0) nhit = -nhit, inside = true;
    let bias: f32 = 1e-4;
    let mut inside: bool = false;
    if raydir.dot(&nhit) > 0.0 {
        nhit = -nhit;
        inside = true;
    }

    //     if ((sphere->transparency > 0 || sphere->reflection > 0) && depth < MAX_RAY_DEPTH) {
    //         float facingratio = -raydir.dot(nhit);
    //         float fresneleffect = mix(pow(1 - facingratio, 3), 1, 0.1);
    //         Vec3f refldir = raydir - nhit * 2 * raydir.dot(nhit);
    //         refldir.normalize();
    //         Vec3f reflection = trace(phit + nhit * bias, refldir, spheres, depth + 1);
    //         Vec3f refraction = 0;
    if sphere.transparency > 0.0 || sphere.reflection > 0.0 && depth < MAX_RAY_DEPTH {
        let facingratio = -raydir.dot(&nhit);
        let fresneleffect = mix((1.0 - facingratio).powf(3.0), 1.0, 0.1);
        let mut refldir = raydir - nhit * Vec3::build1(2.0) * Vec3::build1(raydir.dot(&nhit));
        refldir.normalize();
        let reflection = trace(
            phit + nhit * Vec3::build1(bias),
            refldir,
            spheres,
            depth + 1,
        );
        let mut refraction = Vec3::build1(0.0);
        // was sphere.transparency, nullable?
        //         if (sphere->transparency) {
        //             float ior = 1.1, eta = (inside) ? ior : 1 / ior; // are we inside or outside the surface?
        //             float cosi = -nhit.dot(raydir);
        //             float k = 1 - eta * eta * (1 - cosi * cosi);
        //             Vec3f refrdir = raydir * eta + nhit * (eta *  cosi - sqrt(k));
        //             refrdir.normalize();
        //             refraction = trace(phit - nhit * bias, refrdir, spheres, depth + 1);
        //         }
        if sphere.transparency != 0.0 {
            let ior = 1.1;
            // inside or outside
            let eta = if inside { ior } else { 1.0 / ior };
            let cosi = -nhit.dot(&raydir);
            let k = 1.0 - eta * eta * (1.0 - cosi * cosi);
            let mut refrdir =
                raydir * Vec3::build1(eta) + nhit * Vec3::build1(eta * cosi - k.sqrt());
            refrdir.normalize();
            refraction = trace(
                phit - nhit * Vec3::build1(bias),
                refrdir,
                spheres,
                depth + 1,
            );
        }
        //         // the result is a mix of reflection and refraction (if the sphere is transparent)
        //         surfaceColor = (
        //             reflection * fresneleffect +
        //             refraction * (1 - fresneleffect) * sphere->transparency) * sphere->surfaceColor;
        //     }
        let surface_color = (reflection * Vec3::build1(fresneleffect)
            + refraction * Vec3::build1(1.0 - fresneleffect) * Vec3::build1(sphere.transparency))
            * sphere.surface_color;
    } else {
    //         // it's a diffuse object, no need to raytrace any further
    //         for (unsigned i = 0; i < spheres.size(); ++i) {
    //             if (spheres[i].emissionColor.x > 0) {
    //                 // this is a light
    //                 Vec3f transmission = 1;
    //                 Vec3f lightDirection = spheres[i].center - phit;
    //                 lightDirection.normalize();
    //                 for (unsigned j = 0; j < spheres.size(); ++j) {
    //                     if (i != j) {
    //                         float t0, t1;
    //                         if (spheres[j].intersect(phit + nhit * bias, lightDirection, t0, t1)) {
    //                             transmission = 0;
    //                             break;
    //                         }
    //                     }
    //                 }
    //                 surfaceColor += sphere->surfaceColor * transmission *
    //                 std::max(float(0), nhit.dot(lightDirection)) * spheres[i].emissionColor;
    //             }
    //         }
    //     }
        for vsphere in spheres.iter() {
            
        }
    }

    //     return surfaceColor + sphere->emissionColor;
    // }
    return surface_color + sphere.emission_color;
    // Vec3::build0()
}

fn main() {}
