cmake_minimum_required(VERSION 3.23)
project(homework)

set(CMAKE_CXX_STANDARD 20)

add_executable(homework main.cpp maths_functions/triangulus_Delone/triangulus.cpp maths_functions/triangulus_Delone/triangulus.h include/include.h data_structures/string/GS_string.cpp data_structures/string/GS_string.h data_structures/vector/GS_vector.cpp data_structures/vector/GS_vector.h data_structures/vector/include.h data_structures/string/include.h maths_functions/triangulus_Delone/include.h data_structures/matrix/GS_matrix.cpp data_structures/matrix/GS_matrix.h data_structures/matrix/include.h)
