#include "triangle.h"
#include "swap.h"
#include "display.h"

void draw_filled_triangle(
    int x0, int y0, float z0, float w0,
    int x1, int y1, float z1, float w1,
    int x2, int y2, float z2, float w2,
    uint32_t color)
{

    // sort vert by y-coordinate (y0 < y1 < y2)
    if (y0 > y1)
    {
        int_swap(&y0, &y1);
        int_swap(&x0, &x1);
        float_swap(&z0, &z1);
        float_swap(&w0, &w1);
    }

    if (y1 > y2)
    {
        int_swap(&y1, &y2);
        int_swap(&x1, &x2);
        float_swap(&z1, &z2);
        float_swap(&w1, &w2);
    }

    if (y0 > y1)
    {
        int_swap(&y0, &y1);
        int_swap(&x0, &x1);
        float_swap(&z0, &z1);
        float_swap(&w0, &w1);
    }

    vec4_t point_a = {x0, y0, z0, w0};
    vec4_t point_b = {x1, y1, z1, w1};
    vec4_t point_c = {x2, y2, z2, w2};

    // render upper part flat bottom triangle

    float inv_slope_1 = 0;
    float inv_slope_2 = 0;

    if (y1 - y0 != 0)
    {
        inv_slope_1 = (float)(x1 - x0) / abs(y1 - y0);
    }

    if (y2 - y0 != 0)
    {
        inv_slope_2 = (float)(x2 - x0) / abs(y2 - y0);
    }

    if (y1 - y0 != 0)
    {
        for (int y = y0; y <= y1; y++)
        {
            int x_start = x1 + (y - y1) * inv_slope_1;
            int x_end = x0 + (y - y0) * inv_slope_2;

            if (x_end < x_start)
            {
                int_swap(&x_start, &x_end); // always should be from left to right
            }

            for (int x = x_start; x < x_end; x++)
            {
                // draw_texel(x, y, texture, point_a, point_b, point_c, a_uv, b_uv, c_uv);
                draw_triangle_pixel(x, y, color, point_a, point_b, point_c);
            }
        }
    }

    // render bottom part - flat top triangle

    inv_slope_1 = 0;
    inv_slope_2 = 0;

    if (y2 - y1 != 0)
    {
        inv_slope_1 = (float)(x2 - x1) / abs(y2 - y1);
    }

    if (y2 - y0 != 0)
    {
        inv_slope_2 = (float)(x2 - x0) / abs(y2 - y0);
    }

    if (y2 - y1 != 0)
    {
        for (int y = y1; y <= y2; y++)
        {
            int x_start = x1 + (y - y1) * inv_slope_1;
            int x_end = x0 + (y - y0) * inv_slope_2;

            if (x_end < x_start)
            {
                int_swap(&x_start, &x_end); // always should be from left to right
            }

            for (int x = x_start; x < x_end; x++)
            {
                // draw_texel(x, y, texture, point_a, point_b, point_c, a_uv, b_uv, c_uv);
                draw_triangle_pixel(x, y, color, point_a, point_b, point_c);
            }
        }
    }
}

vec3_t barycentric_weights(vec2_t a, vec2_t b, vec2_t c, vec2_t p)
{
    vec2_t ab = vec2_sub(b, a);
    vec2_t bc = vec2_sub(c, b);
    vec2_t ac = vec2_sub(c, a);
    vec2_t ap = vec2_sub(p, a);
    vec2_t bp = vec2_sub(p, b);

    // calculate the area of a full triangle ABC using cross product of parrallelogram
    float area_triangle_abc = (ab.x * ac.y - ab.y * ac.x);

    // weight alpha is the area of subtriangle BCP divided by area of triangle ABC
    float alpha = (bc.x * bp.y - bp.x * bc.y) / area_triangle_abc;

    // weight beta is the area of subtriangle ACP divided by area of triangle ABC
    float beta = (ap.x * ac.y - ac.x * ap.y) / area_triangle_abc;

    // weight gamma is found by subtracting alpha and beta from 1
    float gamma = 1 - alpha - beta;

    vec3_t weights = {alpha, beta, gamma};
    return weights;
}

void draw_triangle_pixel(
    int x, int y, uint32_t color,
    vec4_t point_a, vec4_t point_b, vec4_t point_c)
{
    vec2_t p = {x, y};
    vec2_t a = vec2_from_vec4(point_a);
    vec2_t b = vec2_from_vec4(point_b);
    vec2_t c = vec2_from_vec4(point_c);

    vec3_t weights = barycentric_weights(a, b, c, p);

    float alpha = weights.x;
    float beta = weights.y;
    float gamma = weights.z;

    // also interpolate the value of 1/w for current pixel
    float interpolated_reciprocal_w = (1 / point_a.w) * alpha + (1 / point_b.w) * beta + (1 / point_c.w) * gamma;

    // adjust 1/w so the pixels that are closer to the camera have smaller values
    interpolated_reciprocal_w = 1.0 - interpolated_reciprocal_w;

    // only draw pixel if value is less than previously stored in z-buffer
    if (interpolated_reciprocal_w < z_buffer[(window_width * y) + x])
    {
        draw_pixel(x, y, color);

        // update the z-buffer with the interpolated value of 1/w
        z_buffer[(window_width * y) + x] = interpolated_reciprocal_w;
    }
};

void draw_texel(
    int x, int y, uint32_t *texture,
    vec4_t point_a, vec4_t point_b, vec4_t point_c,
    tex2_t a_uv, tex2_t b_uv, tex2_t c_uv)
{
    vec2_t p = {x, y};
    vec2_t a = vec2_from_vec4(point_a);
    vec2_t b = vec2_from_vec4(point_b);
    vec2_t c = vec2_from_vec4(point_c);

    vec3_t weights = barycentric_weights(a, b, c, p);

    float alpha = weights.x;
    float beta = weights.y;
    float gamma = weights.z;

    // variables to store the interpolated values of u and v and also 1/w for current pixel
    float interpolated_u;
    float interpolated_v;
    float interpolated_reciprocal_w;

    // perform the interpolation of a;; U/w and V/w using barycentric weights and a factor of 1/w
    interpolated_u = (a_uv.u / point_a.w) * alpha + (b_uv.u / point_b.w) * beta + (c_uv.u / point_c.w) * gamma;
    interpolated_v = (a_uv.v / point_a.w) * alpha + (b_uv.v / point_b.w) * beta + (c_uv.v / point_c.w) * gamma;

    // also interpolate the value of 1/w for current pixel
    interpolated_reciprocal_w = (1 / point_a.w) * alpha + (1 / point_b.w) * beta + (1 / point_c.w) * gamma;

    // divide back both values by 1/w to get the actual values of u and v
    interpolated_u /= interpolated_reciprocal_w;
    interpolated_v /= interpolated_reciprocal_w;

    // map the interpolated values to the texture width and height
    int tex_x = abs((int)(interpolated_u * texture_width)) % texture_width;
    int tex_y = abs((int)(interpolated_v * texture_height)) % texture_height;

    // adjust 1/w so the pixels that are closer to the camera have smaller values
    interpolated_reciprocal_w = 1.0 - interpolated_reciprocal_w;

    // only draw pixel if value is less than previously stored in z-buffer
    if (interpolated_reciprocal_w < z_buffer[(window_width * y) + x])
    {
        draw_pixel(x, y, texture[(texture_width * tex_y) + tex_x]);

        // update the z-buffer with the interpolated value of 1/w
        z_buffer[(window_width * y) + x] = interpolated_reciprocal_w;
    }
}

void draw_textured_triangle(
    int x0, int y0, float z0, float w0, float u0, float v0,
    int x1, int y1, float z1, float w1, float u1, float v1,
    int x2, int y2, float z2, float w2, float u2, float v2,
    uint32_t *texture)
{
    // sort vert by y-coordinate (y0 < y1 < y2)
    if (y0 > y1)
    {
        int_swap(&y0, &y1);
        int_swap(&x0, &x1);
        float_swap(&z0, &z1);
        float_swap(&w0, &w1);
        float_swap(&u0, &u1);
        float_swap(&v0, &v1);
    }

    if (y1 > y2)
    {
        int_swap(&y1, &y2);
        int_swap(&x1, &x2);
        float_swap(&z1, &z2);
        float_swap(&w1, &w2);
        float_swap(&u1, &u2);
        float_swap(&v1, &v2);
    }

    if (y0 > y1)
    {
        int_swap(&y0, &y1);
        int_swap(&x0, &x1);
        float_swap(&z0, &z1);
        float_swap(&w0, &w1);
        float_swap(&u0, &u1);
        float_swap(&v0, &v1);
    }

    // flip the v component to account for inverted UC=coordinates
    v0 = 1.0 - v0;
    v1 = 1.0 - v1;
    v2 = 1.0 - v2;

    // create vector points and tex coordinates after sort

    vec4_t point_a = {x0, y0, z0, w0};
    vec4_t point_b = {x1, y1, z1, w1};
    vec4_t point_c = {x2, y2, z2, w2};
    tex2_t a_uv = {u0, v0};
    tex2_t b_uv = {u1, v1};
    tex2_t c_uv = {u2, v2};

    // render upper part flat bottom triangle

    float inv_slope_1 = 0;
    float inv_slope_2 = 0;

    if (y1 - y0 != 0)
    {
        inv_slope_1 = (float)(x1 - x0) / abs(y1 - y0);
    }

    if (y2 - y0 != 0)
    {
        inv_slope_2 = (float)(x2 - x0) / abs(y2 - y0);
    }

    if (y1 - y0 != 0)
    {
        for (int y = y0; y <= y1; y++)
        {
            int x_start = x1 + (y - y1) * inv_slope_1;
            int x_end = x0 + (y - y0) * inv_slope_2;

            if (x_end < x_start)
            {
                int_swap(&x_start, &x_end); // always should be from left to right
            }

            for (int x = x_start; x < x_end; x++)
            {
                draw_texel(x, y, texture, point_a, point_b, point_c, a_uv, b_uv, c_uv);
            }
        }
    }

    // render bottom part - flat top triangle

    inv_slope_1 = 0;
    inv_slope_2 = 0;

    if (y2 - y1 != 0)
    {
        inv_slope_1 = (float)(x2 - x1) / abs(y2 - y1);
    }

    if (y2 - y0 != 0)
    {
        inv_slope_2 = (float)(x2 - x0) / abs(y2 - y0);
    }

    if (y2 - y1 != 0)
    {
        for (int y = y1; y <= y2; y++)
        {
            int x_start = x1 + (y - y1) * inv_slope_1;
            int x_end = x0 + (y - y0) * inv_slope_2;

            if (x_end < x_start)
            {
                int_swap(&x_start, &x_end); // always should be from left to right
            }

            for (int x = x_start; x < x_end; x++)
            {
                draw_texel(x, y, texture, point_a, point_b, point_c, a_uv, b_uv, c_uv);
            }
        }
    }
}