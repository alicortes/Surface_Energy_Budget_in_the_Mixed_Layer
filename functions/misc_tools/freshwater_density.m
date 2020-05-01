function density = freshwater_density(t,S)
%Function to compute freshwater density using Chen and Millero
%Nature, vol 266 page 708;  S = total dissolved salt in parts per thousand
%density = density at sea level in kg/m3%
%
S = S* .9951437 ;
%
%compute density _______________________

density = 0.999835 + (6.7914e-5.*t) - (9.0894e-6.*t.^2) + (1.0171e-7.*t.^3) ...
	  - ( 1.2846e-9   .*  t.^4) ...
	  + ( 1.1592e-11  .*  t.^5) ...
	  - ( 5.0125e-14  .*  t.^6) ...
	  + S .* (8.221e-4 - ( 3.87e-6  .*  t) ...
		+ ( 4.99e-8  .* t.^2 ));
density=density.*1000;


