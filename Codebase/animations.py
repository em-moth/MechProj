from manim import *

#manim -pql animations.py CreateCircle

class CreateCircle(Scene):
    def construct(self):
        ax = Axes()

        self.play(Create(ax))
        self.wait(3)

class ThreeDExample(MovingCameraScene):
    def construct(self):
        axes = ThreeDAxes()
        
        
