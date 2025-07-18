#include "opengl-framework/opengl-framework.hpp"
#include "utils.hpp"

#include <glm/gtc/constants.hpp>

/*struct Circle
{
    
};*/

struct Parallelogram
{
    glm::vec2 origin;
    glm::vec2 vectorA;
    glm::vec2 vectorB;

    glm::vec2 get_random_point_inside() const {
        float u = utils::rand(0.0f, 1.0f);
        float v = utils::rand(0.0f, 1.0f);
        return origin + u * vectorA + v * vectorB;
    }
};

Parallelogram g_parallelogram;

auto bezier3(const glm::vec2& p0, const glm::vec2& p1, const glm::vec2& p2, const glm::vec2& p3) {
    return [p0, p1, p2, p3](float t) {
        glm::vec2 q0 = (1-t)*p0 + t*p1;
        glm::vec2 q1 = (1-t)*p1 + t*p2;
        glm::vec2 q2 = (1-t)*p2 + t*p3;
        
        glm::vec2 r0 = (1-t)*q0 + t*q1;
        glm::vec2 r1 = (1-t)*q1 + t*q2;
        
        return (1-t)*r0 + t*r1;
    };
}

float distance_squared_to_bezier(float t, glm::vec2 point) {
    glm::vec2 bez = bezier3({-.3f, -.3f}, {-0.2f, -.3f}, gl::mouse_position(), {.8f, -.3f})(t);
    return glm::dot(bez - point, bez - point);
}

float find_closest_t(glm::vec2 point, int iterations = 20, float lr = 0.01f) {
    float t = 0.5f;
    for (int i = 0; i < iterations; ++i) {
        float epsilon = 0.001f;
        float d1 = distance_squared_to_bezier(t - epsilon, point);
        float d2 = distance_squared_to_bezier(t + epsilon, point);
        float grad = (d2 - d1) / (2 * epsilon);
        t -= lr * grad;
        t = glm::clamp(t, 0.f, 1.f);
    }
    return t;
}

glm::vec2 get_bezier_normal(float t) {
    glm::vec2 p0 = {-.3f, -.3f};
    glm::vec2 p1 = {-0.2f, 0.5f};
    glm::vec2 p2 = gl::mouse_position();
    glm::vec2 p3 = {.8f, .5f};

    glm::vec2 d = -3.f * (1 - t) * (1 - t) * p0
                  + 3.f * (1 - t) * (1 - t) * p1
                  - 6.f * (1 - t) * t * p1
                  + 6.f * (1 - t) * t * p2
                  - 3.f * t * t * p2
                  + 3.f * t * t * p3;

    glm::vec2 tangent = glm::normalize(d);
    return glm::vec2(-tangent.y, tangent.x);
}

struct Particle
{
    glm::vec2 position;
    glm::vec2 velocity;

    float mass;
    glm::vec2 acceleration;

    float age;
    float lifetime;
    float initial_radius;

    glm::vec4 startColor;
    glm::vec4 endColor;

    float t;

    Particle(float t_, const std::function<glm::vec2(float)>& curve) : t(t_)
    {
        position.x = utils::rand(-gl::window_aspect_ratio(), gl::window_aspect_ratio());
        position.y = utils::rand(0.8f, 1.f);
        
        //position.x = utils::rand(-0.2f, 0.2f);
        //position.y = utils::rand(-0.5f, 0.5f);

        //position = g_parallelogram.get_random_point_inside();

        //float t = utils::rand(0.0f, 1.0f);
        //position = curve(t);

        float epsilon = 0.001f;
        glm::vec2 tangent = glm::normalize(curve(t + epsilon) - curve(t));
        glm::vec2 normal(-tangent.y, tangent.x);

        if (utils::rand(0.f, 1.f) > 0.5f)
            normal = -normal;

        float speed = utils::rand(0.1f, 0.3f);
        velocity = normal * speed;

        //float angle = utils::rand(0.f, 2.f * glm::pi<float>());
        //float speed = utils::rand(0.1f, 0.3f);

        //velocity.x = cos(angle) * speed;
        //velocity.y = sin(angle) * speed;

        mass = 1.f;
        acceleration = glm::vec2(0.0f);

        age = 0.0f;
        lifetime = utils::rand(20.0f, 50.0f);
        initial_radius = 0.01f;//0.02f * (lifetime * 0.1);

        startColor = glm::vec4(
        utils::rand(0.5f, 1.0f),
        utils::rand(0.5f, 1.0f),
        utils::rand(0.5f, 1.0f),
        1.0f
        );

        endColor = glm::vec4(
        utils::rand(0.5f, 1.0f),
        utils::rand(0.5f, 1.0f),
        utils::rand(0.5f, 1.0f),
        1.0f
        );

    }

    void applyForce(const glm::vec2& force)
    {
        acceleration += force / mass;
    }

    void update(float dt)
    {
        age += dt;
        
        glm::vec2 gravity(0.0f, -0.5f * mass);
        acceleration = glm::vec2(0.0f);
        applyForce(gravity);

        float t_closest = find_closest_t(position);
        glm::vec2 closest_point = bezier3({-.3f, -.3f}, {-0.2f, -.3f}, gl::mouse_position(), {.8f, -.3f})(t_closest);

        float dist = glm::length(position - closest_point);

        glm::vec2 normal = get_bezier_normal(t_closest);

        float strength = glm::clamp(0.5f / (dist * dist + 0.05f), 0.0f, 5.0f);
        glm::vec2 force = normal * strength * 0.5f;

        applyForce(force);

        

        /*glm::vec2 airFriction = -0.5f * velocity;
        //applyForce(airFriction);

        glm::vec2 mousePos = gl::mouse_position();
        glm::vec2 springForce = (mousePos - position) * 0.5f;
        //applyForce(springForce);

        glm::vec2 dir = position - mousePos;
        glm::vec2 vortexForce = glm::vec2(-dir.y, dir.x) * 0.3f;
        //applyForce(vortexForce);*/
        
        velocity += acceleration * dt;
        position += velocity * dt;

        /*glm::vec2 prev_pos = position;

        //float t = utils::rand(0.0f, 1.0f);
        t += 0.1f * dt;
        if (t > 1.0f) t -= 1.0f;

        position = bezier3({-.3f, -.3f}, {-0.2f, 0.5f}, gl::mouse_position(), {.8f, .5f})(t);*/
    }

    float getCurrentRadius() const
    {
        //float fadeDuration = 2.0f;
        //float timeLeft = lifetime - age;
        //float lifeRatio = 1.0f - (age / lifetime);
        //float lifeRatio = glm::clamp(timeLeft / fadeDuration, 0.0f, 1.0f);
        //return initial_radius * lifeRatio;

        float t = glm::clamp((lifetime - age) / 2.f, 0.f, 1.f);
        float bounce = t * (1.f + 0.5f * sin(10.f * glm::pi<float>() * t));
        return initial_radius * bounce;
    }

    glm::vec4 getCurrentColor() const
    {
        float t = glm::clamp(age / lifetime, 0.0f, 1.0f);
        //return glm::mix(startColor, endColor, t);

        float easedT = 3.f * t * t - 2.f * t * t * t;
        return glm::mix(startColor, endColor, easedT);
    }

    bool isDead() const
    {
        return age >= lifetime;
    }

    void bounce(const glm::vec2& collision_point, const glm::vec2& normal, float elasticity = 0.7f) {

        float distance_through_wall = glm::length(position - collision_point);

        velocity = glm::reflect(velocity, normal) * elasticity;

        position = collision_point + glm::normalize(velocity) * distance_through_wall;
    }
};

std::optional<glm::vec2> intersection(glm::vec2 origin1, glm::vec2 end1, glm::vec2 origin2, glm::vec2 end2)
{
    glm::vec2 d1 = end1 - origin1;
    glm::vec2 d2 = end2 - origin2;

    glm::vec2 delta_origin = origin2 - origin1;

    glm::mat2 M(d1, -d2);

    float det = glm::determinant(M);

    const float epsilon = 1e-6f;
    if (std::abs(det) < epsilon) {
        return std::nullopt; 
    }

    glm::vec2 params = glm::inverse(M) * delta_origin;
    float t1 = params.x; 
    float t2 = params.y;

    if (t1 >= 0.0f && t1 <= 1.0f && t2 >= 0.0f && t2 <= 1.0f) {
        return origin1 + t1 * d1; 
    }

    return std::nullopt;
}

void check_and_draw_intersection(glm::vec2 seg1_start, glm::vec2 seg1_end, glm::vec2 seg2_start, glm::vec2 seg2_end)
{
    auto intersection_point = intersection(seg1_start, seg1_end, seg2_start, seg2_end);
    
    if (intersection_point) {
        utils::draw_disk(*intersection_point, 0.05f, glm::vec4(1.0f, 0.0f, 0.0f, 1.0f));
    }
}

glm::vec2 getSegmentNormal(glm::vec2 start, glm::vec2 end)
{
    glm::vec2 direction = glm::normalize(end - start);
    return glm::vec2(-direction.y, direction.x);
}


void draw_parametric(std::function<glm::vec2(float)> const& parametric_func, float t_min = 0.0f, float t_max = 1.0f, int steps = 100)
{
    glm::vec2 prev = parametric_func(t_min);
    for (int i = 1; i <= steps; ++i) {
        float t = t_min + (t_max - t_min) * i / float(steps);
        glm::vec2 curr = parametric_func(t);
        utils::draw_line(prev, curr, 0.005, glm::vec4(1, 1, 1, 1));
        prev = curr;
    }
}

auto test_line = [](float t) {
    return glm::vec2(t, t);
};

auto circle = [](float t) {
    float r = 0.5f;
    return glm::vec2(r * cos(t), r * sin(t));
};

auto heart = [](float t){
    float scale = 0.05;
    float x = scale * 16 *pow(sin(t), 3);
    float y = scale * (13 * cos(t) - 5 * cos(2*t) - 2 * cos(3*t) - cos(4*t));
    return glm::vec2(x, y);
};

auto bezier1(const glm::vec2& p0, const glm::vec2& p1) {
    return [p0, p1](float t) {
        return (1-t)*p0 + t*p1;
    };
}

auto bezier2(const glm::vec2& p0, const glm::vec2& p1, const glm::vec2& p2) {
    return [p0, p1, p2](float t) {
        glm::vec2 q0 = (1-t)*p0 + t*p1;
        glm::vec2 q1 = (1-t)*p1 + t*p2;
        return (1-t)*q0 + t*q1;
    };
}


int binomial(int n, int k) {
    if (k > n) return 0;
    if (k == 0 || k == n) return 1;
    return binomial(n - 1, k - 1) + binomial(n - 1, k);
}

float bernstein(int n, int i, float t) {
    return binomial(n, i) * pow(t, i) * pow(1 - t, n - i);
}

auto bezier1_bernstein(const glm::vec2& p0, const glm::vec2& p1) {
    return [p0, p1](float t) {
        return bernstein(1, 0, t) * p0 + bernstein(1, 1, t) * p1;
    };
}

auto bezier2_bernstein(const glm::vec2& p0, const glm::vec2& p1, const glm::vec2& p2) {
    return [p0, p1, p2](float t) {
        return bernstein(2, 0, t) * p0 
             + bernstein(2, 1, t) * p1 
             + bernstein(2, 2, t) * p2;
    };
}

auto bezier3_bernstein(const glm::vec2& p0, const glm::vec2& p1, const glm::vec2& p2, const glm::vec2& p3) {
    return [p0, p1, p2, p3](float t) {
        return bernstein(3, 0, t) * p0 
             + bernstein(3, 1, t) * p1 
             + bernstein(3, 2, t) * p2 
             + bernstein(3, 3, t) * p3;
    };
}

int main()
{
    gl::init("Particules!");
    gl::maximize_window();
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE);

    g_parallelogram.origin = glm::vec2(-0.4f, 0.0f);
    g_parallelogram.vectorA = glm::vec2(0.6f, -0.4f);
    g_parallelogram.vectorB = glm::vec2(0.8f, 0.6f);

    
    // TODO: create an array of particles
    //std::vector<Particle> particles(100);                    //Create particles

    const int particle_count = 100;
    std::vector<Particle> particles;
    particles.reserve(particle_count);

    auto curve = bezier3({-.3f, -.3f}, {-0.2f, -.3f}, gl::mouse_position(), {.8f, -.3f});

    for (int i = 0; i < particle_count; ++i) {
        float t = i / float(particle_count - 1); // entre 0 et 1 inclus
        particles.emplace_back(t, curve);
    }


    //segment 1
    glm::vec2 seg_start = glm::vec2(-1.f, 0.0f);
    glm::vec2 seg_end   = glm::vec2( 1.f, 0.0f);

    //segment 2
    glm::vec2 seg2_start = glm::vec2(0.0f, -0.5f); 

    float thickness = 0.01f;
    glm::vec4 color = glm::vec4(1, 1, 1, 1);

    while (gl::window_is_open())
    {
        glClearColor(0.f, 0.f, 0.f, 1.f);
        glClear(GL_COLOR_BUFFER_BIT);

        const float dt = gl::delta_time_in_seconds();

        particles.erase(
            std::remove_if(particles.begin(), particles.end(),
                [](const Particle& p) { return p.isDead(); }),
            particles.end()
        );  

        // TODO update particles
        // TODO render particles
        for (auto& particle : particles)
        {
            glm::vec2 prev_pos = particle.position;

            particle.update(dt);

            /*if (auto intersect = intersection(seg_start, seg_end, prev_pos, particle.position)) {
                
                glm::vec2 normal = getSegmentNormal(seg_start, seg_end);

                particle.bounce(*intersect, normal);

                particle.position = *intersect + normal * 0.001f;
            }*/

            utils::draw_disk(
                particle.position,  // Position 
                particle.getCurrentRadius(),  // Size
                //glm::vec4(1.f, 1.f, 1.f, 1.f),   // Color
                particle.getCurrentColor()
            );

            //check_and_draw_intersection(seg_start, seg_end, prev_pos, particle.position);
        };

        //utils::draw_line(seg_start, seg_end, thickness, color);
        //utils::draw_line(seg2_start, gl::mouse_position(), thickness, color);
        
        //check_and_draw_intersection(seg_start, seg_end, seg2_start, gl::mouse_position());
        
        //draw_parametric(test_line, 0.0f, 0.5f, 100);
        //draw_parametric(circle, 0.0f, glm::two_pi<float>(), 100);
        //draw_parametric(heart, 0.0f, glm::two_pi<float>(), 200);

        //draw_parametric(bezier1({-.3f, -.3f}, {0.5f, 0.5f}));
        //draw_parametric(bezier2({-.3f, -.3f}, gl::mouse_position(), {.8f, .5f}));
        draw_parametric(bezier3({-.3f, -.3f}, {-0.2f, -.3f}, gl::mouse_position(), {.8f, -.3f}));

        //draw_parametric(bezier1_bernstein({-.3f, -.3f}, {0.5f, 0.5f}));
        //draw_parametric(bezier2_bernstein({-.3f, -.3f}, gl::mouse_position(), {.8f, .5f}));
        //draw_parametric(bezier3_bernstein({-.3f, -.3f}, {-0.2f, 0.5f}, gl::mouse_position(), {.8f, .5f}));

    }
}