#ifndef SPHERE_H
#define SPHERE_H

#include "hittable.h"
#include "vec3.h"

class sphere : public hittable {
public:
    sphere(point3 _center, double _radius, shared_ptr<material> _material)
            : center(_center), radius(_radius), mat(_material) {}

    bool hit(const ray &r, interval ray_t, hit_record &rec) const override {        /*
         * A sphere is hit if the quadratic equation has real roots.
         * P is just the point on the ray, A + tb
         * A is the origin of the ray, and b is the direction of the ray.
         *      (P-C).(P-C) = (x-Cx)^2 + (y-Cy)^2 + (z-Cz)^2 = r^2
         *      (A + tb - C).(A + tb - C) - r^2 = 0
         *      t^2 * b.b  +  2t * b * (A-C)  +  (A-C).(A-C) - r^2 = 0
         * Which gives the quadratic equation
         */

        //    auto a = dot(r.direction(), r.direction()); // b.b
        //    auto b = 2.0 * dot(oc, r.direction()); // 2b (A-C)
        //    auto c = dot(oc, oc) - radius * radius; // (A-C)^2 - r^2
        //    auto discriminant = b * b - 4 * a * c;
        vec3 oc = r.origin() - center;
        auto a = r.direction().length_squared();
        auto half_b = dot(oc, r.direction());
        auto c = oc.length_squared() - radius * radius;

        auto discriminant = half_b * half_b - a * c;
        if (discriminant < 0) return false;
        auto sqrtd = sqrt(discriminant);

        // Find the nearest root that lies in the acceptable range.
        auto root = (-half_b - sqrtd) / a;
        if (!ray_t.surrounds(root)) {
            root = (-half_b + sqrtd) / a;
            if (!ray_t.surrounds(root))
                return false;
        }

        rec.t = root;
        rec.p = r.at(rec.t);
        vec3 outward_normal = (rec.p - center) / radius;
        rec.set_face_normal(r, outward_normal);
        rec.mat = mat;
//        rec.normal = (rec.p - center) / radius;

        return true;
    }

private:
    point3 center;
    double radius;
    shared_ptr<material> mat;
};

#endif
