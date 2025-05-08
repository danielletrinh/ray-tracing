/**
 * Fundamentals of Computer Graphics
 * @author Danielle Trinh 郑云霄 • 2025403099
 * @date 2025-04-23
 */

#include <iostream>
#include <string>

#include "raytracer.h"

int main(int argc, char* argv[]) {
    try {
        std::string scene_file = "../src/scenes/scene.txt";
        std::string output_file = "../src/scenes/output.bmp";

        Raytracer raytracer;
        raytracer.SetInput(scene_file);
        raytracer.SetOutput(output_file);

        std::cout << "Rendering " << scene_file << "...\n";
        raytracer.Run();
        std::cout << "Output saved to " << output_file << "\n";

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}