Code used for measuring contact angle on x-ray microtomography images, based on the method proposed by Scanziani. et al (2017).

Scanziani, A., Singh, K., Blunt, M. J., & Guadagnini, A. (2017). Automatic method for estimation of in situ effective contact angle from X-ray micro tomography images of two-phase flow in porous media. Journal of Colloid And Interface Science, 496, 51–59. https://doi.org/10.1016/j.jcis.2017.02.005.

This code was adapted from the code provided by Scanziani in the page (https://github.com/alessioscanziani/contact-angle-python) in order to bypass the need for the software Thermofisher Avizo 9.5, used for tracking the contact line and extracting the image perpendicual to the plane. This way the the entire code can be excecuted in python.

The acuracy is tested using benchmark test as shown in the following article:

Zevenbergen, C.I., Bazarin, R.L.M., Siebert, D.N. and Santos, L.O.E., 2023. “A comparative analysis of automated contact angle measurement methods for x-ray microtomographic images of two-phase flow in porous media”. doi: 10.26678/ABCM.COBEM2023.COB2023-1723.
