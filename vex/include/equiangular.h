#include <math.h>

struct rayData{
	float offsets[];
	vector sigma_s[];
	vector sigma_t[];
	vector T[];
	float cdfs[];
	float pdfs[];


	// return values for given offset
	void getValuesFromOffset(float offset, pdf; vector _sigma_s, _T)
		{
		int i=0;
		// TODO binary search
		while (offsets[i+1] < offset  && i<len(offsets) ) i++;
		float k = (offset-offsets[i])/(offsets[i+1]-offsets[i]);
		pdf =  lerp(pdfs[i], pdfs[i+1], k);
		_sigma_s = lerp(sigma_s[i], sigma_s[i+1], k);
		vector _sigma_t = lerp(sigma_t[i], sigma_t[i+1], k);
		_T =  T[i] * exp( - (offset-offsets[i]) * _sigma_t );
		}

	// return values for given random value
	void getValuesFromRandom(float r, pdf, offset; vector _sigma_s, _T)
		{
		int i=0;
		r = lerp(cdfs[0], cdfs[-1], r);
		// TODO binary search
		while (cdfs[i+1] < r) i++;
		float k = (r-cdfs[i])/(cdfs[i+1]-cdfs[i]);
		
		pdf = lerp(pdfs[i], pdfs[i+1], k)/cdfs[-1]/(offsets[i+1]-offsets[i]);
		// TODO exponential interpolation
		offset = lerp(offsets[i], offsets[i+1], k);
		_sigma_s = lerp(sigma_s[i], sigma_s[i+1], k);
		vector _sigma_t = lerp(sigma_t[i], sigma_t[i+1], k);
		_T =  T[i] * exp( - (offset-offsets[i]) * _sigma_t );
		}
}

float misWeight(float a, b){
	// balance heuristic
	//return a/(a+b);

	// power heuristic
	float a2 = a*a;
	return a2/(a2+b*b);
}

float getDensity(vector p;string file){
	vector pos = ptransform("space:current", "space:world", p);
	float density=0;
	for (int i=0; i< nprimitives(file);i++){
		density += volumesample(file,i,pos);
		}
	return density;
	}

//compute Transmittance
vector computeT(vector pos, direction, extinction; float stepsize; string file){

	// TODO pos outside bbox when calculate lit fog
	int result;
	vector min,max,near,far;

	// find BBOX
	getbbox(file, min, max);
	near = ptransform("space:world", pos);
	far = ptransform("space:world", pos + direction);

	// find Intersections
	clip(result, near, far, min, max);
	//near = ptransform( "space:world", "space:current", near);
	far = ptransform( "space:world", "space:current", far);
	float dist = distance(far, pos);
	float offset = rand(SID)*stepsize;
	float d = 0;
	vector dir = normalize(direction);
	vector p = pos+dir*offset;
	while (offset < dist)
	{
		p = pos+dir*offset;
		d += getDensity(p,file);
		offset+=stepsize;
	}
	return exp(-extinction * stepsize * d);
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

// compute light distribution in given point
vector computeLight(vector pos; int lid; vector scattering, extinction; float maxdist, height){

	float h1, h2, h;
	vector attenToPos, attenToLight;
	
	
	h = fit01(pow(height, 0.25), 25, 0.0001);
	h1 = ptransform("space:world", 0).y;
	h2 = ptransform("space:world", pos).y;

	if (abs(h1-h2)<0.001 || h==1){
		attenToPos = exp(-extinction * exp(-h*h1) * length(pos));
	}else{
		attenToPos = exp (-extinction * ( length(pos)*(1-exp(h*(h1-h2))) / (h*exp(h*h1)*(h2-h1))));
		}
		
	vector clr = 0;
	illuminance( pos, {0,0,0}, M_PI, "lightmask", getlightname(lid)){
		if (max(Cl)>0.00001){
			shadow(Cl,pos,L);
			
			h1 = ptransform("space:world", pos+L).y;

			if (abs(h1-h2)<0.001 || h==1){
				attenToLight = exp(-extinction * exp(-h*h1) * length(L));
			}else{
				attenToLight = exp (-extinction * ( length(L)*(1-exp(h*(h1-h2))) / (h*exp(h*h1)*(h2-h1))));
				}
			
			clr += Cl * scattering * exp(-h*h2) * attenToPos * attenToLight *0.5;// right value for uniform phase function
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


void compute(vector outColor,
			outOpacity;
			vector perlight[];
			vector scattering, absorption;
			float maxdist, height, shadowscale;
			int samples;
			vector n_freq, n_off; float n_amp, n_rough; int n_octaves;
			vector equiangular, exponential, debug_weights;
			string file){

	if (file == "")
	{
		vector extinction = absorption + scattering;
		float ext_lum = luminance(extinction);
		float dist = min(maxdist,length(I));
		vector nI = normalize(I);
		vector alldistrib = 0;
		if (height ==1 )
		{
			outOpacity = 1 - exp(-dist * extinction);
		}
		else
		{
			float h1, h2, h;
			h1 = ptransform("space:world",0).y;
			h2 = ptransform("space:world", vector(nI*dist)).y;
			h = fit01(pow(height, 0.25), 25, 0.0001);
			if (abs(h1-h2)<0.001 || h==1){
				outOpacity = 1 - exp(-extinction * exp(-h*h1)*dist);
			}else{
				outOpacity = 1 - exp (-extinction * ( dist*(1-exp(h*(h1-h2))) / (h*exp(h*h1)*(h2-h1))));
				}
		}
			
		// get valid light selection
		string lightmask, categories;
		renderstate( "fog:lightmask", lightmask);
		renderstate( "fog:categories", categories);
		if (categories == "") categories = "*"; 
		
		int lights[]  = getlights("lightmask",lightmask, "categories",categories);
		float weight = 1.0/samples;

		
		// loop for per-light evaluation
		for (int i=0; i < len(lights); i++){

			int lid = lights[i];
			vector pos;
			vector clr;
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
				
				if (height ==1 )
				clr=computeLight(pos, lid, scattering, extinction, maxdist);
				else
				clr=computeLight(pos, lid, scattering, extinction, maxdist, height);
				
				 clr *= simpleNoise(pos, n_freq, n_off, n_amp, n_rough, n_octaves);
				other_pdf = ext_lum/(exp(d1*ext_lum)-exp((d1-dist)*ext_lum));
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
				
				if (height ==1 )
				clr=computeLight(pos, lid, scattering, extinction, maxdist);
				else
				clr=computeLight(pos, lid, scattering, extinction, maxdist, height);
				
				clr *= simpleNoise(pos, n_freq, n_off, n_amp, n_rough, n_octaves);
				other_pdf = D/((thetaB - thetaA)*(D*D+(delta-d2)*(delta-d2)));
				clr /= pdf;
				clr *= (weight * misWeight(pdf, other_pdf));
				exponential += clr;
				distrib += clr;
				debug_weights += {0,1,0} * luminance(clr);

			}
			perlight[i] = distrib;
			alldistrib += distrib;
		}
		debug_weights *= weight;
		debug_weights /= (debug_weights.r + debug_weights.g + debug_weights.b);
		outColor = alldistrib;
	}
	else
	{
	int result;
	vector nI = normalize(I);
	vector min,max,near,far;

	// find BBOX
	getbbox(file, min, max);
	near = ptransform( "space:camera", "space:world", {0,0,0});
	far = ptransform( "space:camera", "space:world", I);

	// find Intersections
	clip(result, near, far, min, max);
	near = ptransform( "space:world", "space:camera", near);
	far = ptransform( "space:world", "space:camera", far);

	// do nothing if ray miss the volume BBOX
	if (!result) return;

	float volumesteprate;
	renderstate("object:volumesteprate",volumesteprate);
	// TODO get right value for each volume primitive
	float diagonal = volumevoxeldiameter(file, 0)*0.577350269;// divide by sqrt(3)
	float stepsize = diagonal/volumesteprate;

	vector extinction = absorption + scattering;
	maxdist = distance(near, far);
	
	rayData data;
	float offset = min(stepsize, maxdist)*rand(SID);
	
	vector transmittance = 1;
	float density = 0;
	float sum_pdfs = 0;
	float pdf;
	
	int additional = distance(far,I)<stepsize;
	
	// create ray dataset for decoupled ray marching
	while ((offset<maxdist) && (max(transmittance) > 0.0001))
	{
		vector pos = near + nI *offset;
		
		density = getDensity(pos, file);

		transmittance *= exp(-density * extinction * stepsize);
		pdf = luminance(transmittance * density * scattering);
		sum_pdfs+=pdf;
		
		push(data.offsets, offset);
		push(data.T, transmittance);
		push(data.sigma_s, density * scattering);
		push(data.sigma_t, density * extinction);
		push(data.pdfs,pdf);
		push(data.cdfs,sum_pdfs);

		offset += stepsize;
	}

	//skip calculation if zero density along ray
	if (len(data.T)<2) return;
	if (transmittance==1) return;
	
	// find min&max distance
	float a,b;
	int i = 0;
	int l = len(data.sigma_s);
	
	while (data.sigma_s[i]==0 && i<l) i++;
	a = data.offsets[i]+length(near);
	i = l-1;
	while (data.sigma_s[i]==0  && i<l) i--;
	b = data.offsets[i]+length(near);

	vector alldistrib = 0;
	vector sigma_s, T;
	
	// get valid light selection
	string lightmask, categories;
	renderstate( "fog:lightmask", lightmask);
	renderstate( "fog:categories", categories);
	if (categories == "") categories = "*"; 

	int lights[]  = getlights("lightmask",lightmask, "categories",categories);
	float weight = 1.0/samples;

	// loop for per-light evaluation
	for (i=0; i < len(lights); i++){

		int lid = lights[i];
		vector pos;
		vector clr;
		vector distrib = 0;
		float other_pdf;
		
		// get Light origin
		vector posL = ptransform( getlightname( lights[i] ) , "space:camera" , {0,0,0});
		
		// get coord of closest point to light along (infinite) ray
		float delta = dot(posL, nI);
		
		// get distance this point is from light
		float D = length(delta*nI - posL);
		
		// get angle of endpoints
		float thetaA = atan( a - delta, D);
		float thetaB = atan( b - delta, D);
		
		for (int sample=0; sample<samples; sample++){
			
			// stratified sampling
			float r = (rand(sample-SID)+sample) / samples;

			// EQUIANGULAR SAMPLPING
			float tt = D*tan(lerp(thetaA, thetaB, r));
			float d1 = delta + tt;
			pdf = D/((thetaB - thetaA)*(D*D + tt*tt));

			pos = nI * d1;
			data->getValuesFromOffset(d1-length(near), other_pdf, sigma_s, T);
			//clr=computeLight(pos, lid, scattering, extinction, maxdist);
			clr=0;
			illuminance(pos, {0,0,0}, "lightmask", getlightname(lid))
			{
				if (max(shadow(Cl, pos, L))>0.0001)
				{
					vector contrib = Cl*T*sigma_s* computeT(pos, L, extinction*shadowscale, stepsize, file)*0.5;
					clr += contrib;
				}
			}

			clr /= pdf;
			clr *= (weight * misWeight(pdf, other_pdf));
			equiangular += clr;
			distrib += clr;
			debug_weights += {1,0,0} * luminance(clr);


			// sampling proportional to sigmaS * delta * Transmittance
			data->getValuesFromRandom(r, pdf, offset, sigma_s, T);
			other_pdf = D/((thetaB - thetaA)*(D*D+(delta-offset-length(near))*(delta-offset-length(near))));
			
			pos = near + nI *offset;

			clr = 0;
			illuminance(pos, {0,0,0}, "lightmask", getlightname(lid))
			{
				if (max(shadow(Cl, pos, L))>0.0001)
				{
					vector contrib = Cl*T*sigma_s* computeT(pos, L, extinction*shadowscale, stepsize, file)*0.5;
					clr += contrib;
				}
			}
			clr/=pdf;
			clr *= (weight * misWeight(pdf, other_pdf));
			exponential += clr;
			distrib += clr;
			debug_weights += {0,1,0} * luminance(clr);

		}
		perlight[i] = distrib;
		alldistrib += distrib;
	}
	debug_weights *= weight;
	debug_weights /= (debug_weights.r + debug_weights.g + debug_weights.b);
	debug_weights.b = len(data.T);
	outColor = alldistrib;


	//if (len(data.T)) ;i
	outOpacity = 1-data.T[-1];
	}
}

