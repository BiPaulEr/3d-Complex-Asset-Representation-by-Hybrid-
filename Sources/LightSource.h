#ifndef LIGHTSOURCE_H
#define LIGHTSOURCE_H

#include <glm/glm.hpp>
#include <glm/ext.hpp>

class LightSource : public Transform {
public:
  LightSource() {}
  LightSource(glm::vec3 position, glm::vec3 color, float intensity,   glm::vec3 direction, float coneAngle_) :
  Transform(), position(position), color(color), intensity(intensity),  direction(direction) ,coneAngle(coneAngle_){}
  ~LightSource() {}
  
  float getIntensity() {return intensity;}
  float getConeAngle() { return coneAngle; }
  void setConeAngle(float coneAngle_) {  coneAngle=coneAngle_; }
  void setPosition(glm::vec3 Position_) { position = Position_; }
  float getAc() {return ac;}
  float getAl() {return al;}
  float getAq() {return aq;}
  void setAc(float ac_) {  ac=ac_; }
  void setAl(float al_) { al=al_; }
  void setAq(float aq_) {  aq=aq_; }
  //permet de détacher la source et la caméra
  void TranformPositionandOrientation(glm::mat4 ViewMatrix) {
	  glm::vec4 pos = ViewMatrix * glm::vec4(position, 1);
	  position = glm::vec3(pos.w);
	  glm::vec4 dir = ViewMatrix * glm::vec4(direction, 1);
	  direction = normalize(glm::vec3(dir.w));
	}

  glm::vec3 getDirection() {return direction;}
  
  void setDirection(glm::vec3 direction_) {direction = direction_;}
  glm::vec3 getPosition() { return position; }
  glm::vec3 getColor() { return color; }
private:
 
	float intensity;
  float ac=0.0f;
  float al= 0.0f;
  glm::vec3 position;
  glm::vec3 color;
  float aq= 0.0f;
  float coneAngle;
  glm::vec3 direction;
};

#endif // LIGHTSOURCE_H
