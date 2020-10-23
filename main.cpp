#include <iostream>
#include <chrono>
#include <thread>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include "simulator.hpp"

int main(){
  Simulator simulator;
  bool success = simulator.init();
  if (!success) return 1;

  volatile bool terminate_thread = false;
  auto simulation_thread = std::thread([&](){
    auto time = std::chrono::system_clock::now();
    while(!terminate_thread && !simulator.isFinished()){
      auto new_time = std::chrono::system_clock::now();
      if (static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(new_time - time).count() / 1000.0) < simulator.get_dt()) continue;
      else time = new_time;
      simulator.step();
    }
  });

  igl::opengl::glfw::Viewer viewer;
  igl::opengl::glfw::imgui::ImGuiMenu menu;
  viewer.plugins.push_back(&menu);
  viewer.core().is_animating = true;
  viewer.core().animation_max_fps = 30.0;
  viewer.data().show_faces = false;
  viewer.data().show_lines = true;
  viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer & viewer) -> bool{
    if(viewer.core().is_animating){
      while(true){
        Eigen::MatrixXd V = simulator.V();
        Eigen::MatrixXi F = simulator.F().cast<int>();
        if (F.maxCoeff()<V.rows()){
          viewer.data().clear();
          viewer.data().set_mesh(V, F);
          break;
        }
      }
    }
    return false;
  };
  viewer.data().set_mesh(simulator.V(), simulator.F().cast<int>());
  viewer.launch();

  terminate_thread = true;
  simulation_thread.join();

  return 0;
}