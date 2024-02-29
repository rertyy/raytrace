#ifndef CAMERA_H
#define CAMERA_H

#include "sphere.h"
#include "hittable_list.h"
#include <fstream>
#include "vec3.h"
#include "ray.h"
#include <iostream>
#include "rtweekend.h"

#include "color.h"
#include "hittable.h"
#include "material.h"

#include "rtweekend.h"

#include "color.h"
#include "hittable.h"

#include <iostream>

class camera {
public:
    double aspect_ratio = 1.0;  // Ratio of image width over height
    int image_width = 100;  // Rendered image width in pixel count
    int samples_per_pixel = 10;   // Count of random samples for each pixel
    int max_depth = 10;   // Maximum number of ray bounces into scene

    double vfov = 90;  // Vertical view angle (field of view)

    void render(const hittable &world) {
        initialize();

        std::clog << image_width << "x" << image_height << " image\n";

        // Use file instead of cout to write to image.ppm directly
        std::ofstream file("image.ppm");

        file << "P3\n" << image_width << ' ' << image_height << "\n255\n";

        for (int j = 0; j < image_height; ++j) {
            std::clog << "\rScanlines remaining: " << (image_height - j) << ' ' << std::flush;
            for (int i = 0; i < image_width; ++i) {

                // Enable anti-aliasing by sampling multiple rays per pixel.
                color pixel_color(0, 0, 0);
                for (int sample = 0; sample < samples_per_pixel; ++sample) {
                    ray r = get_ray(i, j);
                    pixel_color += ray_color(r, max_depth, world);
                }
                write_color(file, pixel_color, samples_per_pixel);

            }
        }
        file.flush();
        std::clog << "\rDone.                 \n";
    }

private:
    int image_height;   // Rendered image height
    point3 center;         // Camera center
    point3 pixel00_loc;    // Location of pixel 0, 0
    vec3 pixel_delta_u;  // Offset to pixel to the right
    vec3 pixel_delta_v;  // Offset to pixel below

    void initialize() {
        image_height = static_cast<int>(image_width / aspect_ratio);
        image_height = (image_height < 1) ? 1 : image_height;

        center = point3(0, 0, 0);

        // Determine viewport dimensions.
        auto focal_length = 1.0;
        auto viewport_height = 2.0;
        auto viewport_width = viewport_height * (static_cast<double>(image_width) / image_height);

        // Calculate the vectors across the horizontal and down the vertical viewport edges.
        auto viewport_u = vec3(viewport_width, 0, 0);
        auto viewport_v = vec3(0, -viewport_height, 0);

        // Calculate the horizontal and vertical delta vectors from pixel to pixel.
        pixel_delta_u = viewport_u / image_width;
        pixel_delta_v = viewport_v / image_height;

        // Calculate the location of the upper left pixel.
        auto viewport_upper_left =
                center - vec3(0, 0, focal_length) - viewport_u / 2 - viewport_v / 2;
        pixel00_loc = viewport_upper_left + 0.5 * (pixel_delta_u + pixel_delta_v);
    }


    color ray_color(const ray &r, int depth, const hittable &world) const {
        hit_record rec;

        // If we've exceeded the ray bounce limit, no more light is gathered.
        if (depth <= 0)
            return {0, 0, 0};

        // Use polymorphism to determine which hit function to call.
        if (world.hit(r, interval(0.001, infinity), rec)) {
            // If interval min=0, precision errors where rays which bounce hit the same object surface again
            // causing self-shadowing. This is called shadow acne. Without using more complicated methods,
            // setting min=0.001 causes the sphere to be lighter, with reduced self-shadowing

            // If the ray hits something in the world,
            // use the normal vector to the sphere to color the pixel.
            // We will get a nice colorful sphere on a green world.
            // The green world is because using the normal vector.
            // (notice the small sphere also has green at the top).
            //  return 0.5 * (rec.normal + color(1, 1, 1));


            // Now using diffuse reflections or refractions instead to get grey spheres.
            // Grey because 0.5 is reflected
            // And recursively do so
            //  vec3 direction = random_on_hemisphere(rec.normal);

            // Now use Lambertian reflection instead, so shadows are more accurately modelled.
            // The underneath of the sphere has darker shadows
            // And the spheres are tinted blue from the sky
//            vec3 direction = rec.normal + random_unit_vector();
//            // You can vary the amount of reflectance by changing the 0.5
//            return 0.7 * ray_color(ray(rec.p, direction), depth - 1, world);

            // Now use the material to determine the color of the pixel.
            ray scattered;
            color attenuation;
            if (rec.mat->scatter(r, rec, attenuation, scattered))
                return attenuation * ray_color(scattered, depth - 1, world);
            return {0, 0, 0};
        }


        // Not used anymore due to bespoke hit function
        //    // If the ray hits the sphere, use the normal vector to the sphere to color the pixel.
        //    if (t > 0.0) {
        //        vec3 N = unit_vector(r.at(t) - centre);
        //        // Normal vectors are between -1 and 1, so map them to 0 and 1.
        //        return 0.5 * color(N.x() + 1, N.y() + 1, N.z() + 1);
        //    }

        // If the ray doesn't hit the sphere, use the sky gradient.
        vec3 unit_direction = unit_vector(r.direction());
        // In a unit vector, each component is between -1 and 1
        // So add 1 to each component and divide by 2 to get a value between 0 and 1.
        auto a = 0.5 * (unit_direction.y() + 1.0);
        // Linearly blend white and blue depending on the y-coordinate of the ray.
        return (1.0 - a) * color(1.0, 1.0, 1.0) + a * color(0.5, 0.7, 1.0);
    }

    ray get_ray(int i, int j) const {
        // Get a randomly sampled camera ray for the pixel at location i,j.

        auto pixel_center = pixel00_loc + (i * pixel_delta_u) + (j * pixel_delta_v);
        auto pixel_sample = pixel_center + pixel_sample_square();

        auto ray_origin = center;
        auto ray_direction = pixel_sample - ray_origin;

        return ray(ray_origin, ray_direction);
    }

    vec3 pixel_sample_square() const {
        // Returns a random point in the square surrounding a pixel at the origin.
        auto px = -0.5 + random_double();
        auto py = -0.5 + random_double();
        return (px * pixel_delta_u) + (py * pixel_delta_v);
    }

    vec3 pixel_sample_disk(double radius) const {
        // Generate a sample from the disk of given radius around a pixel at the origin.
        auto p = radius * random_in_unit_disk();
        return (p[0] * pixel_delta_u) + (p[1] * pixel_delta_v);
    }
};

#endif