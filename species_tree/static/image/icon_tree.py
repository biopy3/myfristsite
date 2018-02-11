import turtle

def branch(length, level):
    if level <= 0:
        return

    turtle.forward(length)
    turtle.left(45)
    branch(0.6 * length, level - 1)
    turtle.right(90)
    branch(0.6 * length, level - 1)
    turtle.left(45)
    turtle.backward(length)
    return

turtle.left(90)
branch(100, 6)
ts = turtle.getscreen()
ts.getcanvas().postscript(file="icon_tree.ps")