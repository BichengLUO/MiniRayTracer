#include <Windows.h>
#include <vector>
#include <opencv2\opencv.hpp>
#define WIDTH 500
#define HEIGHT 500
#define UNIT_PER_PIXEL 0.001
#define PI 3.14159265
#define EPS 0.0001
#define REFLECTION_COUNT 2

struct Color { unsigned char r, g, b; };
struct Vec3d { double x, y, z; };
struct Sphere { Vec3d pos; double r; Color c; double shineness; double reflection; }; //Information about the color, position and size of a sphere
struct Ray { Vec3d pos, dir; }; //A ray shooting from 'pos' with direction vector 'dir'
const double ambient = 0.2; //Ambient color intensity
const char *title = "Mini Ray Tracer (by Bicheng LUO, Tsinghua University)"; // Title of the window

double dot(Vec3d a, Vec3d b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vec3d reflect(Vec3d in, Vec3d normal) {
	double in_nl2 = 2 * dot(in, normal);
	return Vec3d{ in.x - in_nl2 * normal.x, in.y - in_nl2 * normal.y, in.z - in_nl2 * normal.z };
}

Vec3d normalize(Vec3d in) {
	double len = sqrt(dot(in, in));
	return Vec3d{ in.x / len, in.y / len, in.z / len };
}

Vec3d negative(Vec3d in) {
	return Vec3d{-in.x, -in.y, -in.z};
}

//Calculate the intersection between a ray and a sphere
bool ray_intersect(Ray ray, Sphere sphere, Vec3d *its, Vec3d *normal, double *depth) {
	double a = dot(ray.dir, ray.dir);
	double b = 2 * (ray.dir.x * (ray.pos.x - sphere.pos.x) + ray.dir.y * (ray.pos.y - sphere.pos.y) + ray.dir.z * (ray.pos.z - sphere.pos.z));
	double c = (ray.pos.x - sphere.pos.x) * (ray.pos.x - sphere.pos.x) + (ray.pos.y - sphere.pos.y) * (ray.pos.y - sphere.pos.y) + (ray.pos.z - sphere.pos.z) * (ray.pos.z - sphere.pos.z) - sphere.r * sphere.r;
	double delta = b*b - 4 * a * c;
	if (delta < 0) return false; //The ray does not intersect the sphere
	double t1 = (-b + sqrt(delta)) / (2 * a);
	double t2 = (-b - sqrt(delta)) / (2 * a);
	double t;
	if ((t1 < t2 && t1 > EPS) || (t2 < 0 && t1 > EPS)) t = t1;
	else if ((t2 < t1 && t2 > EPS) || (t1 < 0 && t2 > EPS)) t = t2;
	else return false; //The intersection point does not lie on the ray
	its->x = ray.pos.x + ray.dir.x * t;
	its->y = ray.pos.y + ray.dir.y * t;
	its->z = ray.pos.z + ray.dir.z * t;
	normal->x = (its->x - sphere.pos.x) / sphere.r;
	normal->y = (its->y - sphere.pos.y) / sphere.r;
	normal->z = (its->z - sphere.pos.z) / sphere.r;
	*depth = t * sqrt(a);
	return true;
}

bool reflect_color(const std::vector<Sphere*> &scene, const Vec3d &light, const Ray &ray,
	Color *reflect_color, Ray *reflection_ray, Sphere **hit) {
	Vec3d its_r, normal_r; //Intersection point and its normal
	double depth_r = std::numeric_limits<double>::max();
	Sphere *hit_r = nullptr;
	bool is_shadow = false;
	Vec3d its, normal;
	double depth;
	for (Sphere *sphere : scene) {
		if (ray_intersect(ray, *sphere, &its, &normal, &depth) && depth < depth_r) {
			its_r = its;
			normal_r = normal;
			depth_r = depth;
			hit_r = sphere; //Update the nearest hit object with the smaller depth
		}
	}
	if (hit_r != nullptr) {
		*hit = hit_r;
		Ray light_ray{ its_r, Vec3d{ light.x - its_r.x, light.y - its_r.y, light.z - its_r.z } };
		for (Sphere *sphere : scene) {
			if (ray_intersect(light_ray, *sphere, &its, &normal, &depth)) {
				is_shadow = true; //Determine whether the pixel is in the shadow
				break;
			}
		}

		double diffuse_r = MAX(dot(normalize(light_ray.dir), normal_r), 0); //Calculate the diffuse color intensity
		double specular_r = pow(MAX(dot(normalize(negative(reflect(light_ray.dir, normal_r))), normalize(negative(ray.dir))), 0), 20 * hit_r->shineness);
		double intensity_r = (is_shadow ? diffuse_r / 3.0 : diffuse_r + specular_r) + ambient; //Decrease the color intensity for the pixels in shadow
		*reflect_color = Color{ MIN(hit_r->c.b * intensity_r, 255), MIN(hit_r->c.g * intensity_r, 255), MIN(hit_r->c.r * intensity_r, 255) };
		*reflection_ray = Ray{ its_r, reflect(ray.dir, normal_r) };
		return true;
	}
	return false;
}

//The entry point of the program
int WINAPI WinMain(HINSTANCE hinst, HINSTANCE, LPSTR, int) {
	cv::Mat canvas(HEIGHT, WIDTH, CV_8UC3, cv::Scalar(0, 0, 0));
	Vec3d light{10, 10, 10};
	std::vector<Sphere*> scene;
	scene.push_back(new Sphere{ Vec3d{ 0, 0, 0 }, 0.5, Color{ 50, 125, 50 }, 0.3, 6 }); //Add a green sphere
	scene.push_back(new Sphere{ Vec3d{ 0, 0.8, 0.8 }, 0.2, Color{ 255, 100, 100 }, 0.3, 0.1 }); //Add a red and small sphere
	scene.push_back(new Sphere{ Vec3d{ 0.8, 0.8, 0 }, 0.3, Color{ 100, 100, 255 }, 0.3, 0.1 }); //Add a blue and small sphere
	int frame = 0;
	while (true) { //The main loop
		canvas.setTo(cv::Scalar(0, 0, 0)); //Clear the canvas
		scene[1]->pos.x = scene[2]->pos.z = 0.8 * sin(frame / 180.0 * PI);
		scene[1]->pos.y = scene[2]->pos.x = 0.8 * cos(frame / 180.0 * PI);
		scene[1]->pos.z = scene[2]->pos.y = cos((frame++) / 180.0 * PI); //Rotate the small spheres
		Vec3d its, normal;
		double depth;
#pragma omp parellel for
		for (int j = 0; j < HEIGHT; j++) {
			uchar *ptr = canvas.ptr<uchar>(j);
			for (int i = 0; i < WIDTH; i++) {
				double pix_y = (i - WIDTH / 2.0) * UNIT_PER_PIXEL;
				double pix_x = (j - HEIGHT / 2.0) * UNIT_PER_PIXEL;

				Ray ray{ Vec3d{ 0, 0, 5 }, Vec3d{ pix_x, pix_y, -1 } }; //Cast each pixel to the specific ray
				
				Color hit_color[REFLECTION_COUNT];
				double hit_color_weight[REFLECTION_COUNT];
				double total_hit_color_weight = 0;
				memset(hit_color_weight, 0, sizeof(double)*REFLECTION_COUNT);
				Sphere *hit[REFLECTION_COUNT];
				for (int k = 0; k < REFLECTION_COUNT; k++) {
					if (!reflect_color(scene, light, ray, &hit_color[k], &ray, &hit[k]))
						break;
					hit_color_weight[k] = k == 0 ? 1 : hit_color_weight[k - 1] * hit[k - 1]->reflection;
					total_hit_color_weight += hit_color_weight[k];
				}
				for (int k = 0; k < REFLECTION_COUNT; k++) {
					ptr[3 * i] += hit_color[k].b * hit_color_weight[k] / total_hit_color_weight;
					ptr[3 * i + 1] += hit_color[k].g * hit_color_weight[k] / total_hit_color_weight;
					ptr[3 * i + 2] += hit_color[k].r * hit_color_weight[k] / total_hit_color_weight;
				}
			}
		}
		cv::imshow(title, canvas); //Show the canvas in the window
		if (cv::waitKey(1) == 27) break; //Press 'ESC' to quit
	}
}