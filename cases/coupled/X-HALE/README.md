X-HALE model for SHARPy
========================

This model is the one used in the 2019 IFASD paper
"Nonlinear Response of a Very Flexible Aircraft Under Lateral Gust" by
A. del Carre, P. Teixeira, R. Palacios and C. Cesnik (IFASD-2019-090).

As it is right now, running the script "generate_xhale.py" with python
will generate the case to run with `sharpy xhale_ifasd_o.sharpy`. The default
gust is a lateral 1-cos with gradient 7.5 chords and intensity 15% of the
free stream speed.

Annotations in the generation script are provided in the comments.
