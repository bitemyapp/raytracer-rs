extern crate byteorder;

// #if defined __linux__ || defined __APPLE__
// // "Compiled for Linux
// #else
// // Windows doesn't define these values by default, Linux does
// #define M_PI 3.141592653589793
// #define INFINITY 1e8
// #endif
// #define MAX_RAY_DEPTH 5

// use byteorder::{ReadBytesExt, WriteBytesExt, BigEndian, LittleEndian};
// use std::cmp::max;
use std::f32::NAN;
use std::cmp::Ordering;
use std::fs::File;
use std::io::prelude::*;
use std::ops::Add;
use std::ops::Mul;
use std::ops::Neg;
use std::ops::Sub;

const PI: f32 = 3.141592653589793;
const INFINITY: f32 = 1e8;
const MAX_RAY_DEPTH: u64 = 5;

#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
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
            x: self.x.neg(),
            y: self.y.neg(),
            z: self.z.neg(),
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

#[derive(Clone, Copy, PartialEq)]
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
    fn intersect(&self, rayorig: Vec3, raydir: Vec3, mut t0: &mut f32, mut t1: &mut f32) -> bool {
        let l: Vec3 = self.center - rayorig;
        let tca: f32 = l.dot(&raydir);
        if tca < 0.0 {
            return false;
        }
        let d2: f32 = l.dot(&l) - tca * tca;
        if d2 > self.radius2 { return false; }
        let thc: f32 = (self.radius2 - d2).sqrt();
        // let tc_sub = (tca - thc).clone();
        // let tc_add = (tca + thc).clone();
        *t0 = tca - thc;
        *t1 = tca + thc;
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

fn partial_max<T>(x: T, y: T) -> T where T: PartialOrd {
    match x.partial_cmp(&y) {
        None => x,
        Some(Ordering::Less) => y,
        Some(Ordering::Equal) => x,
        Some(Ordering::Greater) => x,
    }
}
// Vec3f trace(
//     const Vec3f &rayorig,
//     const Vec3f &raydir,
//     const std::vector<Sphere> &spheres,
//     const int &depth)

fn trace(rayorig: Vec3, raydir: Vec3, spheres: &Vec<Sphere>, depth: u64) -> Vec3 {
    // {
    //     //if (raydir.length() != 1) std::cerr << "Error " << raydir << std::endl;
    // if raydir.length() != 1.0 {
    //     panic!("raydir length wasn't 1");
    // }
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
    // println!("Survived osphere");
    //     Vec3f surfaceColor = 0; // color of the ray/surfaceof the object intersected by the ray
    //     Vec3f phit = rayorig + raydir * tnear; // point of intersection
    //     Vec3f nhit = phit - sphere->center; // normal at the intersection point
    //     nhit.normalize(); // normalize normal direction
    let mut surface_color = Vec3::build1(0.0);
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
        surface_color =
            (reflection
             * Vec3::build1(fresneleffect)
             + refraction
             * Vec3::build1(1.0 - fresneleffect)
             * Vec3::build1(sphere.transparency))
            * sphere.surface_color;
    } else {
        //         // it's a diffuse object, no need to raytrace any further
        //         for (unsigned i = 0; i < spheres.size(); ++i) {
        //             if (spheres[i].emissionColor.x > 0) {
        for isphere in spheres.iter() {
            //                 // this is a light
            //                 Vec3f transmission = 1;
            //                 Vec3f lightDirection = spheres[i].center - phit;
            //                 lightDirection.normalize();
            let mut transmission = 1.0;
            let mut light_direction = isphere.center - phit;
            light_direction.normalize();
            //                 for (unsigned j = 0; j < spheres.size(); ++j) {
            for jsphere in spheres.iter() {
                //                     if (i != j) {
                //                         float t0, t1;
                //                         if (spheres[j].intersect(phit + nhit * bias, lightDirection, t0, t1)) {
                //                             transmission = 0;
                //                             break;
                //                         }
                //                     }
                if isphere != jsphere {
                    let mut t0 = NAN;
                    let mut t1 = NAN;
                    if jsphere.intersect(phit + nhit * Vec3::build1(bias),
                                         light_direction,
                                         &mut t0,
                                         &mut t1) {
                        transmission = 0.0;
                        break;
                    }
                }
            }
            // surfaceColor += sphere->surfaceColor * transmission *
            // std::max(float(0), nhit.dot(lightDirection)) * spheres[i].emissionColor;
            let new_surface_color =
                sphere.surface_color * Vec3::build1(transmission)
                * partial_max(Vec3::build1(0.0),
                              Vec3::build1(nhit.dot(&light_direction)))
                * isphere.emission_color;

            surface_color = surface_color + new_surface_color;
        }
    //             }
    //         }
    //     }
    }

    //     return surfaceColor + sphere->emissionColor;
    // }
    // println!("{:?}", surface_color + sphere.emission_color);
    return surface_color + sphere.emission_color;
    // Vec3::build0()
}

// //[comment]
// // Main rendering function. We compute a camera ray for each pixel of the image
// // trace it and return a color. If the ray hits a sphere, we return the color of the
// // sphere at the intersection point, else we return the background color.
// //[/comment]
// void render(const std::vector<Sphere> &spheres)
// {
fn render(spheres: Vec<Sphere>) -> std::io::Result<()> {
//     unsigned width = 640, height = 480;
    const width: usize = 640;
    const height: usize = 480;
//     Vec3f *image = new Vec3f[width * height], *pixel = image;
//     float invWidth = 1 / float(width), invHeight = 1 / float(height);
//     float fov = 30, aspectratio = width / float(height);
//     float angle = tan(M_PI * 0.5 * fov / 180.);
    // let mut image: [Vec3; width * height];
    let mut image: Vec<Vec3> = Vec::new();
    let inv_width = 1.0 / width as f32;
    let inv_height = 1.0 / height as f32;
    let fov = 30.0;
    let aspect_ratio = width as f32 / height as f32;
    let angle = (PI * 0.5 * fov / 180.0).tan();
    let mut pixel = image;
//     // Trace rays
//     for (unsigned y = 0; y < height; ++y) {
//         for (unsigned x = 0; x < width; ++x, ++pixel) {
//             float xx = (2 * ((x + 0.5) * invWidth) - 1) * angle * aspectratio;
//             float yy = (1 - 2 * ((y + 0.5) * invHeight)) * angle;
//             Vec3f raydir(xx, yy, -1);
//             raydir.normalize();
//             *pixel = trace(Vec3f(0), raydir, spheres, 0);
//         }
//     }
    for y in 0..height {
        for x in 0..width {
            let xx = (2.0 * ((x as f32 + 0.5) * inv_width) - 1.0) * angle * aspect_ratio;
            let yy = (1.0 - 2.0 * ((y as f32 + 0.5) * inv_height)) * angle;
            let mut raydir = Vec3 { x: xx, y: yy, z: -1.0 };
            raydir.normalize();
            // pixel[x] = trace(Vec3::build1(0.0), raydir, &spheres, 0);
            pixel.push(trace(Vec3::build1(0.0), raydir, &spheres, 0));
        }
    }
    //     // Save result to a PPM image (keep these flags if you compile under Windows)
    //     std::ofstream ofs("./untitled.ppm", std::ios::out | std::ios::binary);
    //     ofs << "P6\n" << width << " " << height << "\n255\n";
    //     for (unsigned i = 0; i < width * height; ++i) {
    //         ofs << (unsigned char)(std::min(float(1), image[i].x) * 255) <<
    //                (unsigned char)(std::min(float(1), image[i].y) * 255) <<
    //                (unsigned char)(std::min(float(1), image[i].z) * 255);
    //     }
    //     ofs.close();
    //     delete [] image;
    // }
    let mut wtr = vec![];
    // wtr.write_u16::<LittleEndian>(517).unwrap();
    // wtr.write_u16::<LittleEndian>(768).unwrap();
    // wtr.write_u64::<LittleEndian>(width as u64).unwrap();
    // wtr.write(width);
    // wtr.write_u64::<LittleEndian>(height as u64).unwrap();
    wtr.write(b"P6\n");
    write!(wtr, "{}", width);
    wtr.write(b" ");
    write!(wtr, "{}", height);
    wtr.write(b"\n255\n");
    // buffer.write(b"P6\n");
    // buffer.write(width);
    // buffer.write(b" ");
    // buffer.write(height);
    // buffer.write(b"\n255\n");
    // write!(&mut buffer, "{:b}", width);
    // buffer.write(height);
    // write!(&mut buffer, "{:b}", height);
    for i in 0..(width * height) {
        let one: f32 = 1.0;
        // let write_x = one.min(image[i].x) * 255.0;
        // let write_y = one.min(image[i].y) * 255.0;
        // let write_z = one.min(image[i].z) * 255.0;
        // {
        //     println!("{:?} {:?}", pixel[i].x, one.min(pixel[i].x));
        // }
        let write_x = one.min(pixel[i].x) * 255.0;
        let write_y = one.min(pixel[i].y) * 255.0;
        let write_z = one.min(pixel[i].z) * 255.0;
        wtr.write(
            &[write_x as u8,
              write_y as u8,
              write_z as u8,
            ]
        );
        // wtr.write(write_y as u8);
        // wtr.write(write_z as u8);        
        // wtr.write_f32::<LittleEndian>(write_x).unwrap();
        // wtr.write_f32::<LittleEndian>(write_y).unwrap();
        // wtr.write_f32::<LittleEndian>(write_z).unwrap();
        // buffer.write(one.min(image[i].x) * 255.0);
        // buffer.write(one.min(image[i].y) * 255.0);
        // buffer.write(one.min(image[i].z) * 255.0);
        // write!(&mut buffer, "{:b}", one.min(image[i].x) * 255.0);
        // write!(&mut buffer, "{:b}", one.min(image[i].y) * 255.0);
        // write!(&mut buffer, "{:b}", one.min(image[i].z) * 255.0);
    }
    // println!("{:?}", wtr.len());
    let mut buffer = File::create("./untitled.ppm")?;
    buffer.write(&wtr);
    Ok(())
}

// //[comment]
// // In the main function, we will create the scene which is composed of 5 spheres
// // and 1 light (which is also a sphere). Then, once the scene description is complete
// // we render that scene, by calling the render() function.
// //[/comment]
// int main(int argc, char **argv)
// {
//     srand48(13);
//     std::vector<Sphere> spheres;
fn main() {
    let mut spheres = Vec::new();
    spheres.push(
        // build(c: Vec3, r: f32, sc: Vec3, ec: Vec3, refl: f32, transp: f32)
        //     // position, radius, surface color, reflectivity, transparency, emission color
        //     spheres.push_back(Sphere(Vec3f( 0.0, -10004, -20), 10000, Vec3f(0.20, 0.20, 0.20), 0, 0.0));
        Sphere::build(
            Vec3 { x: 0.0, y: -10004.0, z: -20.0 },
            10000.0,
            Vec3 { x: 0.20, y: 0.20, z: 0.20 },
            Vec3::build1(0.0),
            0.0,
            0.0,
        )
    );

    spheres.push(
        //     spheres.push_back(Sphere(Vec3f( 0.0,      0, -20),     4, Vec3f(1.00, 0.32, 0.36), 1, 0.5));
        Sphere::build(
            Vec3 { x: 0.0, y: 0.0, z: -20.0 },
            4.0,
            Vec3 { x: 1.00, y: 0.32, z: 0.36 },
            Vec3::build1(1.0),
            0.5,
            0.0,
        )
    );

    spheres.push(
        //     spheres.push_back(Sphere(Vec3f( 5.0,     -1, -15),     2, Vec3f(0.90, 0.76, 0.46), 1, 0.0));
        Sphere::build(
            Vec3 { x: 5.0, y: -1.0, z: -15.0 },
            2.0,
            Vec3 { x: 0.90, y: 0.76, z: 0.46 },
            Vec3::build1(1.0),
            0.0,
            0.0,
        )
    );

    spheres.push(
        //     spheres.push_back(Sphere(Vec3f( 5.0,      0, -25),     3, Vec3f(0.65, 0.77, 0.97), 1, 0.0));
        Sphere::build(
            Vec3 { x: 5.0, y: 0.0, z: -25.0 },
            3.0,
            Vec3 { x: 0.65, y: 0.77, z: 0.97 },
            Vec3::build1(1.0),
            0.0,
            0.0,
        )
    );

    spheres.push(
        //     spheres.push_back(Sphere(Vec3f(-5.5,      0, -15),     3, Vec3f(0.90, 0.90, 0.90), 1, 0.0));
        Sphere::build(
            Vec3 { x: -5.5, y: 0.0, z: -15.0 },
            3.0,
            Vec3 { x: 0.90, y: 0.90, z: 0.90 },
            Vec3::build1(1.0),
            0.0,
            0.0,
        )
    );


    spheres.push(
        //     // light
        //     spheres.push_back(Sphere(Vec3f( 0.0,     20, -30),     3, Vec3f(0.00, 0.00, 0.00), 0, 0.0, Vec3f(3)));
        Sphere::build(
            Vec3 { x: 0.0, y: 20.0, z: -30.0 },
            3.0,
            Vec3 { x: 0.00, y: 0.00, z: 0.00 },
            Vec3::build1(1.0),
            0.0,
            // Vec3::build1(3.0),
            3.0,
        )
    );
//     render(spheres);
    render(spheres);
//     return 0;
// }
}
