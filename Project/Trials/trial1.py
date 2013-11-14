import cv
from PIL import Image
import numpy

im = cv.LoadImage("sample.png",1)
cv.NamedWindow("example",cv.CV_WINDOW_AUTOSIZE)
cv.ShowImage("example",im)
cv.WaitKey(10000)
cv.DestroyWindow("example")

