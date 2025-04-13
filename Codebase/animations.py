from manim import *


class CreateCircle(Scene):
    def construct(self):
        circle = Circle()  # create a circle
        circle.set_fill(PINK, opacity=0.5)  # set the color and transparency
        self.play(Create(circle))  # show the circle on screen
        t = Text("Hello").to_edge(UL, buff=0.5)
        self.play(Create(t))

        self.wait()
