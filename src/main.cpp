#include "opengl-framework/opengl-framework.hpp"
#include "utils.hpp"

#include <glm/gtc/constants.hpp>

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


    Particle()
    {
        position.x = utils::rand(-gl::window_aspect_ratio(), gl::window_aspect_ratio());
        position.y = utils::rand(-1.f, 1.f);

        float angle = utils::rand(0.f, 2.f * glm::pi<float>());
        float speed = utils::rand(0.1f, 0.3f);

        velocity.x = cos(angle) * speed;
        velocity.y = sin(angle) * speed;

        mass = 1.f;
        acceleration = glm::vec2(0.0f);

        age = 0.0f;
        lifetime = utils::rand(20.0f, 50.0f);
        initial_radius = 0.02f * (lifetime * 0.1);

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
        //applyForce(gravity);

        glm::vec2 airFriction = -0.5f * velocity;
        //applyForce(airFriction);

        glm::vec2 mousePos = gl::mouse_position();
        glm::vec2 springForce = (mousePos - position) * 0.5f;
        //applyForce(springForce);

        glm::vec2 dir = position - mousePos;
        glm::vec2 vortexForce = glm::vec2(-dir.y, dir.x) * 0.3f;
        //applyForce(vortexForce);
        
        velocity += acceleration * dt;
        position += velocity * dt;
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

glm::vec2 getSegmentNormal(glm::vec2 start, glm::vec2 end) {
    glm::vec2 direction = glm::normalize(end - start);
    return glm::vec2(-direction.y, direction.x);
}

int main()
{
    gl::init("Particules!");
    gl::maximize_window();
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE);

    // TODO: create an array of particles
    std::vector<Particle> particles(100);

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

            if (auto intersect = intersection(seg_start, seg_end, prev_pos, particle.position)) {
                
                glm::vec2 normal = getSegmentNormal(seg_start, seg_end);

                particle.bounce(*intersect, normal);

                particle.position = *intersect + normal * 0.001f;
            }

            utils::draw_disk(
                particle.position,  // Position 
                particle.getCurrentRadius(),  // Size
                //glm::vec4(1.f, 1.f, 1.f, 1.f),   // Color
                particle.getCurrentColor()
            );

            check_and_draw_intersection(seg_start, seg_end, prev_pos, particle.position);
        };

        utils::draw_line(seg_start, seg_end, thickness, color);
        utils::draw_line(seg2_start, gl::mouse_position(), thickness, color);
        
        check_and_draw_intersection(seg_start, seg_end, seg2_start, gl::mouse_position());
        
        
    }
}