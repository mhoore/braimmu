#version 3.6;
global_settings {assumed_gamma 1.0}
#default{finish{ambient 0.2 diffuse 0.6 specular 0.2 roughness .1}}
#include "colors.inc"
#include "textures.inc"

//camera {location <40.0, 60.0, 80.0>
//	look_at  <40.0, 60.0, 0.0>
//	up    <0,1,0>
//	right  <-1.,0,0>
//	sky <0,1,0>}

//camera{ orthographic angle  33
//	location <6.5, 40.0, 60.0>
//	look_at  <6.5, 15.0, 0.0>
//   right -x*image_width/image_height}

//light_source { <60,20,100> color White}
//light_source { <20,20,100> color White}

background { rgb <1.0, 1.0, 1.0>}

///colors
#declare mem = rgbft <0.6, 0.6, 0.0, 0.0, 0.0>;

#declare r_bead = 2.0; // *** Note: r_bead = sigma / 2 ***

//Textures
#declare MEM = texture{pigment{color mem}finish {ambient 0.2 diffuse 0.6 specular 0.2 roughness .1}}

//Data
