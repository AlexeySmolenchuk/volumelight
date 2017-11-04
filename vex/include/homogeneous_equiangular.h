#include <math.h>


float misWeight(float a, b){
//    // balance heuristic
//    return a/(a+b);

    // power heuristic
    float a2 = a*a;
    return a2/(a2+b*b);

    // float a3 = a*a*a;
    // return a3/(a3+b*b*b);
}

// compute light distribution in given point
vector computeLight(vector pos; int lid; vector scattering, extinction; float maxdist){
    vector clr = 0;
    illuminance( pos, {0,0,0}, M_PI, "lightmask", getlightname(lid)){
        if (max(Cl)>0.00001){
            shadow(Cl,pos,L);
            
            clr += Cl * scattering * exp(-min(maxdist, (length(L) + length(pos))) * extinction) *0.5;// right value for uniform phase function
        }
    }
    return clr;
}

// this function defined like noise inside LitFog shader
float simpleNoise(vector pos ,n_freq, n_off; float n_amp, n_rough; int n_octaves){
    float nval=0;
    if (n_amp != 0 && n_octaves > 0)
    {
        vector PN = wo_space(pos) * n_freq - n_off;
        float scale = n_amp;
        nval = n_amp * (noise(PN) - 0.5);
        if (n_octaves > 1)
        {
            for (int o = n_octaves; o-- > 1; )
            {
                scale *= n_rough;
                PN *= 2.0;
                // Here, we have 3 octaves of noise
                nval += scale * (noise(PN) - 0.5);
            }
        }
    }
    return 1.0 + max(0, nval);
}


void compute(vector outColor, outOpacity; vector perlight[]; vector scattering, absorption; float maxdist; int samples;
                vector n_freq, n_off; float n_amp, n_rough; int n_octaves;
                vector equiangular,exponential,uniform,debug_weights){
    
    vector extinction = absorption + scattering;
    float ext_lum = luminance(extinction);
    float dist = min(maxdist,length(I));
    vector alldistrib = 0;
    outOpacity = exp(-dist * extinction);
    
    // get valid light selection
    string lightmask, categories;
    renderstate( "fog:lightmask", lightmask);
    renderstate( "fog:categories", categories);
    if (categories == "") categories = "*"; 
    
    int lights[]  = getlights("lightmask",lightmask, "categories",categories);
    float weight = 1.0/samples;
    vector nI = normalize(I);
    
    // loop for per-light evaluation
    for (int i=0; i < len(lights); i++){

        int lid = lights[i];
        vector pos;
        vector clr;
        float scale;
        vector distrib = 0;
        float pdf, other_pdf;
        
        // get Light origin
        vector posL = ptransform( getlightname( lights[i] ) , "space:camera" , {0,0,0});
        
        // get coord of closest point to light along (infinite) ray
        float delta = dot(posL, nI);
        
        // get distance this point is from light
        float D = length(delta*nI - posL);
        
        // get angle of endpoints
        float thetaA = atan(0.0 - delta, D);
        float thetaB = atan(dist - delta, D);
        
        for (int sample=0; sample<samples; sample++){
            
            // stratified sampling
            float r = (rand(sample-SID)+sample) / samples;

            // EQUIANGULAR SAMPLPING
            float tt = D*tan(lerp(thetaA, thetaB, r));
            float d1 = delta + tt;
            pdf = D/((thetaB - thetaA)*(D*D + tt*tt));

            pos = nI * d1;
            clr=computeLight(pos, lid, scattering, extinction, maxdist);
            clr *= simpleNoise(pos, n_freq, n_off, n_amp, n_rough, n_octaves);
            other_pdf = ext_lum/(exp(d1*ext_lum)-exp((d1-dist)*ext_lum));
            // other_pdf = pdf3;
            clr /= pdf;
            clr *= (weight * misWeight(pdf, other_pdf));
            equiangular += clr;
            distrib += clr;
            debug_weights += {1,0,0} * luminance(clr);



            // EXPONENTIAL SAMPLING
            // remap u to account for finite max distance
            float minU = exp(-ext_lum*dist);
            float a1 = r * (1.0 - minU) + minU;

            // sample with pdf proportional to exp(-sig*d)
            float d2 = -log(a1) / ext_lum;
            pdf  = ext_lum * a1 / (1.0 - minU);

            pos = nI * d2;
            clr = computeLight(pos, lid, scattering, extinction, maxdist);
            clr *= simpleNoise(pos, n_freq, n_off, n_amp, n_rough, n_octaves);
            other_pdf = D/((thetaB - thetaA)*(D*D+(delta-d2)*(delta-d2)));
            // other_pdf = 1.0/dist;
            clr /= pdf;
            clr *= (weight * misWeight(pdf, other_pdf));
            exponential += clr;
            distrib += clr;
            debug_weights += {0,1,0} * luminance(clr);



//            // UNIFORM SAMPLING
//            float d3 = dist*r;
//            pdf = 1.0/dist;

//            pos = normalize(I)*d3;
//            clr=computeLight(pos, lid, scattering, extinction);
//            //other_pdf = D/((thetaB - thetaA)*(D*D+(delta-d3)*(delta-d3)));
//            //other_pdf = ext_lum/(exp(d3*ext_lum)-exp((d3-dist)*ext_lum));
//            //clr*=dist*misWeight(pdf3,other_pdf);
//            clr/=pdf;
//            clr *= weight;
//            distrib += clr;
//            //uniform += clr*weight;
//            //debug_weights+={0,0,1}*misWeight(pdf3, pdf1);
        }
        perlight[i] = distrib;
        alldistrib += distrib;
    }
    debug_weights *= weight;
    debug_weights /= (debug_weights.r + debug_weights.g + debug_weights.b);
    outColor = alldistrib;
}
