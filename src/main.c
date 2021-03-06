#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <SDL2/SDL.h>
#include "upng.h"
#include "array.h"
#include "display.h"
#include "vector.h"
#include "matrix.h"
#include "light.h"
#include "camera.h"
#include "triangle.h"
#include "texture.h"
#include "mesh.h"

#define MAX_TRIANGLES_PER_MESH 10000
triangle_t triangles_to_render[MAX_TRIANGLES_PER_MESH];
int num_triangles_to_render = 0;

mat4_t world_matrix;
mat4_t proj_matrix;
mat4_t view_matrix;

bool is_running = false;
int previous_frame_time = 0;
float delta_time = 0;

void setup(void)
{
    render_method = RENDER_TEXTURED;
    cull_method = CULL_BACKFACE;

    // allocate required mem in bytes to hold the color buffer and z buffer
    color_buffer = (uint32_t *)malloc(sizeof(uint32_t) * window_width * window_height);
    z_buffer = (float *)malloc(sizeof(float) * window_width * window_height);

    // create SDL texture that is used to display color buffer
    color_buffer_texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA32, SDL_TEXTUREACCESS_STREAMING, window_width, window_height);

    // initialize the perspective projection matrix
    float fov = M_PI / 2;
    float aspect = (float)window_height / (float)window_width;
    float znear = 0.1;
    float zfar = 100.0;

    proj_matrix = mat4_make_perspective(fov, aspect, znear, zfar);

    // loads vertex and face values for the mesh data structure
    load_obj_file_data("./assets/f22.obj");
    load_png_texture_data("./assets/f22.png");
}

void process_input(void)
{
    SDL_Event event;
    SDL_PollEvent(&event);

    switch (event.type)
    {
    case SDL_QUIT:
        is_running = false;
        break;

    case SDL_KEYDOWN:
        if (event.key.keysym.sym == SDLK_ESCAPE)
            is_running = false;
        if (event.key.keysym.sym == SDLK_1)
            render_method = RENDER_WIRE_VERTEX;
        if (event.key.keysym.sym == SDLK_2)
            render_method = RENDER_WIRE;
        if (event.key.keysym.sym == SDLK_3)
            render_method = RENDER_FILL_TRIANGLE;
        if (event.key.keysym.sym == SDLK_4)
            render_method = RENDER_FILL_TRIANGLE_WIRE;
        if (event.key.keysym.sym == SDLK_5)
            render_method = RENDER_TEXTURED;
        if (event.key.keysym.sym == SDLK_6)
            render_method = RENDER_TEXTURED_WIRE;
        if (event.key.keysym.sym == SDLK_c)
            cull_method = CULL_BACKFACE;
        if (event.key.keysym.sym == SDLK_x)
            cull_method = CULL_NONE;
        if (event.key.keysym.sym == SDLK_UP)
            camera.position.y += 3.0 * delta_time;
        if (event.key.keysym.sym == SDLK_DOWN)
            camera.position.y -= 3.0 * delta_time;
        if (event.key.keysym.sym == SDLK_a)
            camera.yaw += 1.0 * delta_time;
        if (event.key.keysym.sym == SDLK_d)
            camera.yaw -= 1.0 * delta_time;
        if (event.key.keysym.sym == SDLK_w)
        {
            camera.forward_velocity = vec3_mul(camera.direction, 5.0 * delta_time);
            camera.position = vec3_add(camera.position, camera.forward_velocity);
        }
        if (event.key.keysym.sym == SDLK_s)
        {
            camera.forward_velocity = vec3_mul(camera.direction, 5.0 * delta_time);
            camera.position = vec3_sub(camera.position, camera.forward_velocity);
        }
        break;
    }
}

// update function frame by frame with fixed time step
void update(void)
{
    int time_to_wait = FRAME_TARGET_TIME - (SDL_GetTicks() - previous_frame_time);

    if (time_to_wait > 0 && time_to_wait <= FRAME_TARGET_TIME)
    {
        SDL_Delay(time_to_wait);
    }

    // get delta time converted to seconds
    delta_time = (SDL_GetTicks() - previous_frame_time) / 1000.0;

    previous_frame_time = SDL_GetTicks();

    // initilize the counter of triangles to render for the current frame
    num_triangles_to_render = 0;

    // change the mesh scale, rotation and translation values per animation frame
    mesh.translation.z = 5.0;

    // find the target
    vec3_t target = {0.0, 0.0, 1.0};
    mat4_t camera_yaw_rotation = mat4_make_rotation_y(camera.yaw);
    camera.direction = vec3_from_vec4(mat4_mul_vec4_project(camera_yaw_rotation, vec4_from_vec3(target)));

    target = vec3_add(camera.position, camera.direction);
    vec3_t up_direction = {0.0, 1.0, 0.0};

    view_matrix = mat4_look_at(camera.position, target, up_direction);

    // create scale, rotation and translation matrices
    mat4_t scale_matrix = mat4_make_scale(mesh.scale.x, mesh.scale.y, mesh.scale.z);
    mat4_t translation_matrix = mat4_make_translation(mesh.translation.x, mesh.translation.y, mesh.translation.z);
    mat4_t rotation_matrix_x = mat4_make_rotation_x(mesh.rotation.x);
    mat4_t rotation_matrix_y = mat4_make_rotation_y(mesh.rotation.y);
    mat4_t rotation_matrix_z = mat4_make_rotation_z(mesh.rotation.z);

    // loop all triangles of the mesh
    int num_faces = array_length(mesh.faces);
    for (int i = 0; i < num_faces; i++)
    {
        face_t mesh_face = mesh.faces[i];

        vec3_t face_vertices[3];
        face_vertices[0] = mesh.vertices[mesh_face.a];
        face_vertices[1] = mesh.vertices[mesh_face.b];
        face_vertices[2] = mesh.vertices[mesh_face.c];

        vec4_t transformed_vertices[3];

        // loop all three vertices of current face and apply transforms
        for (int j = 0; j < 3; j++)
        {
            vec4_t transformed_vertex = vec4_from_vec3(face_vertices[j]);

            // update world matrix combining scale, rotation and translation matrices
            world_matrix = mat4_identity();

            // order: first scale, than rotate, than translate
            world_matrix = mat4_t_mul_mat4(scale_matrix, world_matrix);
            world_matrix = mat4_t_mul_mat4(rotation_matrix_z, world_matrix);
            world_matrix = mat4_t_mul_mat4(rotation_matrix_y, world_matrix);
            world_matrix = mat4_t_mul_mat4(rotation_matrix_x, world_matrix);
            world_matrix = mat4_t_mul_mat4(translation_matrix, world_matrix);

            // multiplying the world matrix by original vector
            transformed_vertex = mat4_t_mul_vec4(world_matrix, transformed_vertex);

            // multiply the view matrix by the transformed vertex
            transformed_vertex = mat4_t_mul_vec4(view_matrix, transformed_vertex);

            // save transformed vertex to array of transformed vertices
            transformed_vertices[j] = transformed_vertex;
        }

        // Backface culling to toest if current face should be projected
        // if (cull_method == CULL_BACKFACE)

        vec3_t vector_a = vec3_from_vec4(transformed_vertices[0]); /*    A   */
        vec3_t vector_b = vec3_from_vec4(transformed_vertices[1]); /*   / \  */
        vec3_t vector_c = vec3_from_vec4(transformed_vertices[2]); /*  B---C */

        // get vector substractions of A-B and A-C
        vec3_t vector_ab = vec3_sub(vector_b, vector_a);
        vec3_t vector_ac = vec3_sub(vector_c, vector_a);
        vec3_normalize(&vector_ab);
        vec3_normalize(&vector_ac);

        // compute the face normal (cross product to find perpendicular vector)
        vec3_t normal = vec3_cross(vector_ab, vector_ac);
        vec3_normalize(&normal);

        vec3_t origin = {0, 0, 0};

        // find the vector between a point in the triangle (A) and the camera position
        vec3_t camera_ray = vec3_sub(origin, vector_a);

        // calculate how aligned is the camera ray to the face normal using dot-product
        float dot_normal_camera = vec3_dot(normal, camera_ray);

        if (cull_method == CULL_BACKFACE)
        {
            // bypass triangles that are looking away from the camera
            if (dot_normal_camera < 0)
            {
                continue;
            }
        }

        vec4_t projected_points[3];

        // loop all three vertices of current face and project them to 2D screen coordinates
        for (int j = 0; j < 3; j++)
        {
            // project the current vertex
            projected_points[j] = mat4_mul_vec4_project(proj_matrix, transformed_vertices[j]);

            // scale to the current view
            projected_points[j].x *= (window_width / 2.0);
            projected_points[j].y *= (window_height / 2.0);

            // invert the y values to account for flipped screen y coordinates
            projected_points[j].y *= -1;

            // translate to the center of the screen
            projected_points[j].x += (window_width / 2.0);
            projected_points[j].y += (window_height / 2.0);
        }

        // calculate the shade intensity based on how aligned is the face normal and the light ray direction
        float light_intensity_factor = -vec3_dot(normal, light.direction);

        // calculate triangle color based on light angle
        uint32_t triangle_color = light_apply_intensity(mesh_face.color, light_intensity_factor);

        triangle_t projected_triangle = {
            .points = {
                {projected_points[0].x, projected_points[0].y, projected_points[0].z, projected_points[0].w},
                {projected_points[1].x, projected_points[1].y, projected_points[1].z, projected_points[1].w},
                {projected_points[2].x, projected_points[2].y, projected_points[2].z, projected_points[2].w}},
            .texcoords = {
                {mesh_face.a_uv.u, mesh_face.a_uv.v},
                {mesh_face.b_uv.u, mesh_face.b_uv.v},
                {mesh_face.c_uv.u, mesh_face.c_uv.v},
            },
            .color = triangle_color};

        // save projected triangle to the array of triangles to render
        // array_push(triangles_to_render, projected_triangle);
        if (num_triangles_to_render < MAX_TRIANGLES_PER_MESH)
        {
            triangles_to_render[num_triangles_to_render] = projected_triangle;
            num_triangles_to_render++;
        }
    }
}

void render(void)
{
    draw_grid();

    // loop all projected triangles
    for (int i = 0; i < num_triangles_to_render; i++)
    {
        triangle_t triangle = triangles_to_render[i];

        // draw filled triangle
        if (render_method == RENDER_FILL_TRIANGLE || render_method == RENDER_FILL_TRIANGLE_WIRE)
        {
            draw_filled_triangle(
                triangle.points[0].x, triangle.points[0].y, triangle.points[0].z, triangle.points[0].w,
                triangle.points[1].x, triangle.points[1].y, triangle.points[1].z, triangle.points[1].w,
                triangle.points[2].x, triangle.points[2].y, triangle.points[2].z, triangle.points[2].w,
                triangle.color);
        }

        // draw textured triangle
        if (render_method == RENDER_TEXTURED || render_method == RENDER_TEXTURED_WIRE)
        {
            draw_textured_triangle(
                triangle.points[0].x, triangle.points[0].y, triangle.points[0].z, triangle.points[0].w, triangle.texcoords[0].u, triangle.texcoords[0].v, // Vertex A
                triangle.points[1].x, triangle.points[1].y, triangle.points[1].z, triangle.points[1].w, triangle.texcoords[1].u, triangle.texcoords[1].v, // Vertex B
                triangle.points[2].x, triangle.points[2].y, triangle.points[2].z, triangle.points[2].w, triangle.texcoords[2].u, triangle.texcoords[2].v, // Vertex C
                mesh_texture);
        }

        // draw wireframe triangle
        if (render_method == RENDER_WIRE || render_method == RENDER_WIRE_VERTEX || render_method == RENDER_FILL_TRIANGLE_WIRE || render_method == RENDER_TEXTURED_WIRE)
        {
            draw_triangle(triangle.points[0].x, triangle.points[0].y, triangle.points[1].x, triangle.points[1].y, triangle.points[2].x, triangle.points[2].y, 0xFFFFFFFF);
        }

        if (render_method == RENDER_WIRE_VERTEX)
        {
            draw_rect(triangle.points[0].x - 3, triangle.points[0].y - 3, 6, 6, 0xFFFF0000);
            draw_rect(triangle.points[1].x - 3, triangle.points[1].y - 3, 6, 6, 0xFFFF0000);
            draw_rect(triangle.points[2].x - 3, triangle.points[2].y - 3, 6, 6, 0xFFFF0000);
        }
    }

    render_color_buffer();

    clear_color_buffer(0xFF000000);
    clear_z_buffer();

    SDL_RenderPresent(renderer);
}

void free_resources(void)
{
    free(color_buffer);
    free(z_buffer);
    upng_free(png_texture);
    array_free(mesh.faces);
    array_free(mesh.vertices);
}

int main(void)
{
    is_running = initialize_window();

    setup();

    while (is_running)
    {
        process_input();
        update();
        render();
    }

    destroy_window();
    free_resources();

    return 0;
}