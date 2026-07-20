################################################################################
## IMPORT MODULES
################################################################################
import numpy as np
import pygame as pg

################################################################################
## DEFINE PARAMETERS
################################################################################
# Game parameters
NAME = 'Level Control'
VERSION = 'Dev: Draw'

# internal units
# Distance: m
# Time: s

# unit conversions
k_to_u = 1.e3 # e.g. kg_to_g
m_to_px = 100.

# margins
TOP_MARGIN = 100
BOTTOM_MARGIN = 100
LEFT_MARGIN = 100
RIGHT_MARGIN = 100

# tank parameters
TANK_DIM = np.array([400,400],dtype=int) # pixels
TANK_LEFT = LEFT_MARGIN
TANK_RIGHT = TANK_LEFT + TANK_DIM[0]
TANK_TOP = TOP_MARGIN
TANK_BOTTOM = TANK_TOP + TANK_DIM[1]
TANK_COL = 'yellow2'
TANK_WIDTH = 5

# liquid parameters
h = TANK_DIM[1]
LIQ_COL = 'lightblue'

# screen parameters
SCREEN_DIM = TANK_DIM + np.array([LEFT_MARGIN+RIGHT_MARGIN,TOP_MARGIN+BOTTOM_MARGIN])
SCREEN_COL = 'black'
FPS = 60

# flow parameters


################################################################################
## DEFINE CLASSES
################################################################################


################################################################################
## MAIN CODE
################################################################################
# initialize pygame
pg.init()
clock = pg.time.Clock()
surface = pg.display.set_mode(SCREEN_DIM)
pg.display.set_caption(f'{NAME:s} {VERSION:s}')

# animation loop
running = True
alive = True
while running:

    ############################################################################
    # input
    ############################################################################
    for event in pg.event.get():
        # quit
        if event.type == pg.QUIT or (event.type == pg.KEYUP and event.key == pg.K_q):
            running = False

    ############################################################################
    # dynamics
    ############################################################################
    if h >= 0:
        pass
    else:
        Vout = 0

    if h <= TANK_DIM[1]:
        pass
    else: 
        Vin = 0

    Fin = Vin
    Fout = Vout*0.1*np.sqrt(h)
    h += Fin - Fout

    ############################################################################
    # render the screen
    ############################################################################
    surface.fill(SCREEN_COL)
    # liquid
    pos = (TANK_LEFT, TANK_BOTTOM - h)
    dim = (TANK_DIM[0], h)
    liq = pg.Rect(pos, dim)
    pg.draw.rect(surface, LIQ_COL, liq)


    # tank
    pg.draw.lines(surface, TANK_COL, False,TANK,TANK_WIDTH)
    #outlet
    pos = (LEFT_MARGIN + TANK_DIM[0]//2 - VALVE//2 * Fout, TANK_BOTTOM-TANK_WIDTH//2)
    dim = (VALVE *Fout, BOTTOM_MARGIN + TANK_WIDTH//2)
    out = pg.Rect(pos, dim)
    pg.draw.rect(surface, LIQ_COL, out)


    pg.display.flip()

    clock.tick(FPS)

pg.display.quit()
pg.quit()
