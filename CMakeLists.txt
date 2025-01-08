cmake_minimum_required(VERSION 3.27)
project(MyCppCode_Vscode)

set(CMAKE_CXX_STANDARD 23)

add_executable(MyCppCode_Vscode main.cpp)

target_compile_definitions(MyCppCode_Vscode PRIVATE LOCAL)