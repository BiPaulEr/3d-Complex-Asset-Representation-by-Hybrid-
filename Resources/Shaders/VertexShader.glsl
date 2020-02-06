#version 450 core // Minimal GL version support expected from the GPU

layout(location=0) in vec3 vPosition; // The 1st input attribute is the position (CPU side: glVertexAttrib 0)
layout(location=1) in vec3 vNormal;
layout(location=2) in vec2 vTexCoord;

//Component
layout(location=3) in float vCurrentComponent;

//Deuxieme Pos
layout(location=4) in vec3 vPosition2;

uniform mat4 projectionMat, modelViewMat, normalMat; 
uniform float Mode;
uniform float Scale;
uniform float Interpolate;
out vec3 fPosition;
out vec3 fNormal;
out vec2 fTexCoord;
out vec3 fPosition2;

//Component
out float fCurrentComponent ;

void main() {
	vec3 tmp;
	tmp=vPosition;

	// to made interpolation between vertex if needed
	if (Mode <-2 || Interpolate!=0){
		if (vCurrentComponent<0.f){
		}
		else{
			tmp = (vPosition*(1-Scale)+(vPosition2*Scale)) ;	
		}
	}

	
	vec4 p = modelViewMat * vec4 (tmp, 1.0);
    gl_Position =  projectionMat * p; // mandatory to fire rasterization properly
    vec4 n = normalMat * vec4 (vNormal, 1.0);
    fPosition = p.xyz;
    fNormal = normalize (n.xyz);

	//Component transmission
    fTexCoord = vTexCoord;
	fCurrentComponent=vCurrentComponent;
	fPosition2=vPosition2;

}