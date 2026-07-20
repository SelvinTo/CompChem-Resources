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
VERSION = 'Dev: Control Loop'

# margins
TOP_MARGIN = 100
BOTTOM_MARGIN = 100
LEFT_MARGIN = 100
RIGHT_MARGIN = 0

# tank parameters
TANK_DIM = np.array([400,400],dtype=int) # pixels
TANK_LEFT = LEFT_MARGIN
TANK_RIGHT = TANK_LEFT + TANK_DIM[0]
TANK_TOP = TOP_MARGIN
TANK_BOTTOM = TANK_TOP + TANK_DIM[1]
TANK = [None]*4
TANK[0] = (TANK_LEFT,TANK_TOP)
TANK[1] = (TANK_LEFT,TANK_BOTTOM)
TANK[2] = (TANK_RIGHT,TANK_BOTTOM)
TANK[3] = (TANK_RIGHT,TANK_TOP)
TANK_COL = 'yellow2'
TANK_WIDTH = 5

# liquid parameters
h = 0.
LIQ_COL = 'lightblue'

# screen parameters
SCREEN_DIM = TANK_DIM + np.array([LEFT_MARGIN+RIGHT_MARGIN,TOP_MARGIN+BOTTOM_MARGIN])
SCREEN_COL = 'black'
FPS = 60

# flow parameters
VALVE = 50
Vout = 0.
Vin  = 0.5

# control parameters
hSet = 200.
K = 1.
I = 1.
cI = 0.
ubar = 50.

################################################################################
## DEFINE CLASSES
################################################################################
class Slider(object):
    
    num = 0
    # track parms
    track_dx = 150
    # slider parms
    dx = 100
    dy = 50
    # groove parms
    groove_dy = SCREEN_DIM[1] - dy*3
    groove_y = (dy,dy+groove_dy)
    # slider bounds
    ymax = dy // 2 # top position
    ymin = groove_y[-1] - dy // 2 # bottom position

    def __init__(self,lbl,fmt,units,
                      val,val_min,val_max):
        self.idx = Slider.num
        # label
        self.lbl = lbl
        self.fmt = fmt
        self.units = units
        # groove
        self.groove_x = TANK_RIGHT + self.track_dx // 2 + self.track_dx * self.idx
        # value
        self.val = val
        self.val_min = val_min
        self.val_max = val_max
        # slider
        self.click = False
        self.x = self.groove_x - self.dx//2
        self.y = self.ymin + (self.val - self.val_min)/(self.val_max - self.val_min)*-self.groove_dy
        self.rect = pg.Rect(self.x,self.y,self.dx,self.dy)
        self.ddy = 0 # mouse offset

        Slider.num += 1

    def move(self,pos):
        if self.click:
            px,py = pos
            if py+self.ddy < self.ymax:
                self.y = self.ymax
            elif py+self.ddy > self.ymin:
                self.y = self.ymin
            else:
                self.y = py + self.ddy
            self.rect = pg.Rect(self.x,self.y,self.dx,self.dy)
            self.val = self.val_min + (self.y - self.ymin)/-self.groove_dy * (self.val_max - self.val_min)
        
    def draw(self,surface):
        # groove
        groove = ((self.groove_x,self.groove_y[0]),(self.groove_x,self.groove_y[1]))
        pg.draw.lines(surface,'white',False,groove,5)
        # slider
        if not self.click:
            pg.draw.rect(surface,'gray75',self.rect)
        else:
            pg.draw.rect(surface,'gray60',self.rect)
        # label
        text = font12.render(self.lbl, True, 'white')
        textRect = text.get_rect()
        textRect.center = (self.groove_x,self.ymin + self.dy*1.5)
        surface.blit(text,textRect)
        # value
        text = font12.render(self.fmt % self.val, True, 'white')
        textRect = text.get_rect()
        textRect.midright = (self.groove_x,self.ymin + self.dy*2.0)
        surface.blit(text,textRect)
        # units
        text = font12.render(' ' + self.units, True, 'white')
        textRect = text.get_rect()
        textRect.midleft = (self.groove_x,self.ymin + self.dy*2.0)
        surface.blit(text,textRect)

################################################################################
## MAIN CODE
################################################################################
# initialize elements
sliders = []
sliders.append(Slider('Inlet Valve','%.2f','',Vin,0,1))
sliders.append(Slider('Set Point','%.0f','cm',hSet,0,400))
sliders.append(Slider('Instant Gain','%.2f','',K,0.1,10))    ## Added slider for instant K ##
sliders.append(Slider('Integrated Gain','%.2f','',I,0.1,10)) ## Added slider for integral I ##

# screen dimensions
RIGHT_MARGIN = Slider.num * Slider.track_dx
SCREEN_DIM += np.array([RIGHT_MARGIN,0])

# initialize pygame
pg.init()
clock = pg.time.Clock()
surface = pg.display.set_mode(SCREEN_DIM)
pg.display.set_caption(f'{NAME:s} {VERSION:s}')
font12 = pg.font.Font(None,2*12)

# animation loop
running = True
alive = True
while running:

    ############################################################################
    # input
    ############################################################################
    px,py = pg.mouse.get_pos()
    for event in pg.event.get():
        # quit
        if event.type == pg.QUIT or (event.type == pg.KEYUP and event.key == pg.K_q):
            running = False
        if event.type == pg.MOUSEBUTTONDOWN:
            alive = False
            for slider in sliders:
                if slider.rect.collidepoint((px,py)):
                    slider.click = True
                    slider.ddy = slider.y - py
                    if 'Set Point' in slider.lbl:
                        cI = 0
        if event.type == pg.MOUSEBUTTONUP:
            alive = True
            for slider in sliders:
                slider.click = False

    for slider in sliders:
        slider.move((px,py))
        if 'Inlet Valve' in slider.lbl:
            Vin = slider.val
        elif 'Set Point' in slider.lbl:
            hSet = slider.val

    ############################################################################
    # dynamics
    ############################################################################
    error = hSet - h
    cI += error
    u = -(K*error + I*cI)
    Vout = (u + ubar)/(2*ubar)
    if Vout > 1:
        Vout = 1
    elif Vout < 0:
        Vout = 0

    if h > TANK_DIM[1]:
        sliders[0].click = True
        sliders[0].ddy = 0
        sliders[0].move((0,Slider.ymin))
        sliders[0].click = False

    Fin  = Vin
    Fout = Vout*0.1*np.sqrt(h)
    h += Fin - Fout
    ############################################################################
    # render the screen
    ############################################################################
    surface.fill(SCREEN_COL)
    # liquid in tank
    pos = (TANK_LEFT,TANK_BOTTOM-h)
    dim = (TANK_DIM[0],h)
    liq = pg.Rect(pos,dim)
    pg.draw.rect(surface,LIQ_COL,liq)
    # tank
    pg.draw.lines(surface,TANK_COL,False,TANK,TANK_WIDTH)
    # outlet
    pos = (LEFT_MARGIN + TANK_DIM[0]//2 - VALVE//2 * Fout, TANK_BOTTOM - TANK_WIDTH // 2)
    dim = (VALVE * Fout, BOTTOM_MARGIN + TANK_WIDTH // 2)
    out = pg.Rect(pos,dim)
    pg.draw.rect(surface,LIQ_COL,out)
    # inlet
    pos = (LEFT_MARGIN + TANK_DIM[0]//2 - VALVE//2 * Fin,0)
    dim = (VALVE * Fin, TOP_MARGIN + TANK_DIM[1] - h)
    inp = pg.Rect(pos,dim)
    pg.draw.rect(surface,LIQ_COL,inp)
    # level
    level = ((TANK_LEFT,TANK_BOTTOM-hSet),(TANK_RIGHT,TANK_BOTTOM-hSet))
    pg.draw.lines(surface,'white',False,level,TANK_WIDTH)
    # sliders
    for slider in sliders:
        slider.draw(surface)

    text = font12.render('Liquid Level: %.2f cm' % (h), True, 'white')  ## Display current liquid level ##
    textRect = text.get_rect()
    textRect.topleft = (BOTTOM_MARGIN-50, TANK_LEFT + 440)
    surface.blit(text,textRect)

    pg.display.flip()

    clock.tick(FPS)

pg.display.quit()
pg.quit()
