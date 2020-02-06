#version 450 core // Minimal GL version support expected from the GPU

// ensemble des caractéristiques des sources de Lumière
struct LightSource {
	vec3 position;
	vec3 color;
	float intensity;
	float ac;
	float al;
	float aq;
	float coneAngle;
	vec3 direction;
};
//les textures possibles des matériaux
struct Material {
	sampler2D albedoTex;
	sampler2D roughnessTex;
  sampler2D metallicTex;
	sampler2D ambientTex;
	sampler2D toonTex;
};

uniform Material material;

uniform float Scale;
uniform float Interpolate;

in vec3 fPosition; // Shader input, linearly interpolated by default from the previous stage (here the vertex shader)
//normal du sommet
in vec3 fNormal;
in vec2 fTexCoord;

// Component
in float fCurrentComponent;
in vec3 fPosition2;
out vec4 colorResponse; // Shader output: the color response attached to this fragment


uniform LightSource lightSource1;
//le mode tOOn XTOON ou models à micro-facettes
uniform float Mode;

// besoin pour le XTOON
uniform float z_Min;
uniform float z_Max;



void main() {
//mise en place des vecteurs utiles pour le traitement
	vec3 n = normalize (fNormal); // Linear barycentric interpolation does not preserve unit vectors
	vec3 wo = normalize (-fPosition);
	vec3 radiance = vec3 (0.0, 0.0, 0.0);


	//Monte Carlo 

	if ((Mode == 0.f )|| (Mode<-1.5f)) { 
	vec3 wi = normalize(lightSource1.position - fPosition);
	//angle entre la direction de la lumière et le vecteur de la camera à la fPosition
	float cosAngle = clamp(dot(normalize(-wi), normalize(lightSource1.direction)), 0, 1); 
	//est ce que le sommet est éclairé
	
		vec3 wh = normalize(wi+wo);
		

		//distribution GGX
		float alpha = texture(material.metallicTex, fTexCoord).r; //roughness
		float D = pow(alpha, 2)/(3.14*pow(1+dot(pow(alpha, 2)-1, pow(dot(n, wh), 2)), 2));
		
		

		//terme de Fresnel
		vec3 F_0 = vec3(texture(material.metallicTex, fTexCoord).rgb);
		vec3 F = F_0+(1-F_0)*pow((1-max(0, dot(wi, wh))), 5);
		
		vec3 Li = lightSource1.color * lightSource1.intensity;
		float d = distance(lightSource1.position, fPosition);
		//modélisation de l'atténuation lumineuse modulée
		Li = (1/(lightSource1.ac+lightSource1.al*d+lightSource1.aq*d*d))*Li;

		// GGX Geo
		float G = dot(n, wi)/(dot(n, wi)*(1-alpha*sqrt(2/3.14))+alpha*sqrt(2/3.14))*dot(n, wo)/(dot(n, wo)*(1-alpha*sqrt(2/3.14))+alpha*sqrt(2/3.14));
		//terme spéculaire
		vec3 fs = (D*F*G)/(4*dot(n, wi)*dot(n, wo));

		//terme de diffusion
		float fd = texture(material.metallicTex, fTexCoord).r / (6*3.14);	
		//resultance de la radiance si l'objet est dans le cone d'éclairage
		radiance = max(Li * (fs+fd)*texture(material.roughnessTex, fTexCoord).rgb * dot(n, wi), 0);
		
	}
	else if (Mode == 1.f) { //TOON 
		if (dot(n, wo) < 0.3) { //contour
			radiance = vec3 (0.0, 0.0, 0.0);
		}
		else if (dot(n, wo) > 0.8) { //specular spot
			radiance = vec3 (1.0, 1.0, 1.0);
		}
		else { 
		//le reste
			radiance = vec3 (0.0, 1.0, 0.0);
		}
	}
	else if (Mode == 2.f){ 
	//XTOON
		float Value = 1-log(-fPosition.z/z_Min)/log(z_Max/z_Min);

		radiance = texture(material.toonTex, vec2(dot(n, wo), Value)).rgb;	
	}

	// red or green to show the real or potentiel macro surfaces
	else if ((Mode == 3.f) ||(Mode == 4.f) ){
	float r = (fCurrentComponent);
	if (r>0.5f)  {
	radiance = vec3(0.0f,1.0f,0.0f);
	}
	if (r<0.5f){
	radiance = vec3(1.0f,0.0f,0.0f);
	}
	
	if (r<0.f)  {
	radiance = vec3(0.0f,0.0f,0.0f);
	}

	//use if you want to show component related
	/*while (r>=15){
	r=r-15;
	}
	if ( r ==-1 ){
		radiance = vec3(0.0f,0.0f,1.0f);
		}
		if ( r ==0 ){
		radiance = vec3(1.0f,0.0f,0.0f);
		}
		if ( r ==1 ){
		radiance = vec3(0.0f,0.0f,1.0f);
		}
		if ( r ==2 ){
		radiance = vec3(0.0f,1.0f,0.0f);
		}
		if ( r ==3 ){
		radiance = vec3(1.0f,1.0f,0.0f);
		}
		if ( r ==5 ){
		radiance = vec3(0.0f,1.0f,1.0f);
		}
		if ( r ==6){
		radiance = vec3(1.0f,0.0f,1.0f);
		}
		if ( r ==7 ){
		radiance = vec3(0.5f,0.0f,0.0f);
		}
		if ( r ==8 ){
		radiance = vec3(0.0f,0.0f,0.5f);
		}
		if ( r ==9 ){
		radiance = vec3(0.0f,0.5f,0.0f);
		}
		if ( r ==10 ){
		radiance = vec3(0.5f,1.0f,0.0f);
		}
		if ( r ==11 ){
		radiance = vec3(0.0f,1.0f,0.5f);
		}
		if ( r ==12){
		radiance = vec3(0.5f,0.0f,1.0f);
		}
		if ( r ==13 ){
		radiance = vec3(1.0f,0.5f,0.0f);
		}
		if ( r ==14 ){
		radiance = vec3(0.0f,0.5f,1.0f);
		}
		if ( r ==15){
		radiance = vec3(1.0f,0.0f,0.5f);
		}*/


	}
	// la réponse finale
	if (Mode<-2.f || Interpolate!=0){
		if (fCurrentComponent<-1.2f){
			if(fCurrentComponent<-1.7f){
			colorResponse = vec4 (radiance, Scale); 
			}else{
			colorResponse = vec4 (radiance, 1-Scale); 
			}
			
		}else{
		colorResponse = vec4 (radiance, 1.f);
		}
	}
	else{
	colorResponse = vec4 (radiance, 1.f);
	}
	
}
