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
        lifetime = utils::rand(2.0f, 5.0f);
        initial_radius = 0.02f * lifetime;

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
};

int main()
{
    gl::init("Particules!");
    gl::maximize_window();
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE);

    // TODO: create an array of particles
    std::vector<Particle> particles(100);
    

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
            particle.update(dt);

            utils::draw_disk(
                particle.position,  // Position 
                particle.getCurrentRadius(),  // Size
                //glm::vec4(1.f, 1.f, 1.f, 1.f),   // Color
                particle.getCurrentColor()
            );
        };
    }
}