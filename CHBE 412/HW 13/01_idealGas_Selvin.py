################################################################################
## IMPORT MODULES
################################################################################
import numpy as np
import pygame as pg
rng = np.random.default_rng()

################################################################################
## DEFINE PARAMETERS
################################################################################
# Game parameters
NAME = 'Ideal Gas'
VERSION = 'Dev: Temperature'

# internal units
# Mass: g/mol = atomic mass unit (amu)
# Amount: unit
# Temperature: K
# Distance: nm
# Time: ps
# Energy: amu*nm^2/ps^2
# Pressure: amu/nm/ps^2

# unit conversions
bar_to_Pa = 1.e5
k_to_u = 1.e3 # e.g. kg_to_g
u_to_n = 1.e9 # e.g. m_to_nm
u_to_p = 1.e12
u_to_z = 1.e21
mol_to_u = 6.02214e23
nm_to_px = 20. ## changed resolution to 20 pixels/nm ##
m_to_px = u_to_n * nm_to_px
Jmol_to_Eint = (u_to_n/u_to_p)**2*k_to_u

# universal constants
kB = 8.3144626 * Jmol_to_Eint # energy/atom/K

# system parameters
N = np.array([2,2,2])
T = 300
TIMESTEP = 0.1

# atom parameters
COLORS = ['red','orange','yellow2','green','blue','purple','magenta','cyan']

# box parameters
BOX_DIM = np.array([500,500],dtype=int) # pixels
BOX_POS = np.array([0,0],dtype=int)
BOX_COL = 'black'
box = np.array(BOX_DIM / nm_to_px, dtype=float) # nm

# screen parameters
BOTTOM_MARGIN = 0
RIGHT_MARGIN = 0
SCREEN_DIM = BOX_DIM + np.array([RIGHT_MARGIN,BOTTOM_MARGIN])
SCREEN_COL = 'darkblue'
FPS = 100

################################################################################
## DEFINE CLASSES
################################################################################
class Component(object):

    def __init__(self,mass,diameter,color=None):
        self.mass = mass
        self.diameter = diameter
        self.rad = diameter/2
        if color:
            self.color = color
        else:
            self.color = rng.choice(COLORS)
            COLORS.remove(self.color)

    def init_pos(self,num,box):
        self.num = num
        self.pos = rng.random((self.num,2))*(box-self.diameter) + self.rad

    def init_vel(self,T):
        self.vel = rng.normal(scale=np.sqrt(kB*T/self.mass),size=(self.num,2))

    def move(self,dt,box):
        for iAtom in range(self.num):
            self.pos[iAtom] += self.vel[iAtom]*dt
            for iDim in range(2):
                if ((self.pos[iAtom,iDim] > box[iDim]-self.rad and self.vel[iAtom,iDim] > 0) or 
                    (self.pos[iAtom,iDim] < self.rad and self.vel[iAtom,iDim] < 0)):
                    self.vel[iAtom,iDim] *= -1

    def draw(self,surface):
        rad = round(self.rad*nm_to_px)
        for iAtom in range(self.num):
            pos = np.round(self.pos[iAtom]*nm_to_px)
            pg.draw.circle(surface,self.color,pos,rad)

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
        self.groove_x = BOX_DIM[0] + self.track_dx // 2 + self.track_dx * self.idx
        # value
        self.val = val
        self.val_min = val_min
        self.val_max = val_max
        # slider
        self.click = False
        self.x = self.groove_x - self.dx//2
        self.y = self.ymin + (self.val - self.val_min)/(self.val_max - self.val_min)*-self.groove_dy
        self.rect = pg.Rect(self.x,self.y,self.dx,self.dy)
        self.ddy = 0 # offset from mouse click to top left corner

        Slider.num += 1

    def move(self,pos):
        if slider.click:
            px,py = pos
            if py+self.ddy < self.ymax:
                self.y = self.ymax
            elif py+slider.ddy > self.ymin:
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
        slider = pg.Rect(self.x,self.y,self.dx,self.dy)
        if not self.click:
            pg.draw.rect(surface,'gray75',slider)
        else:
            pg.draw.rect(surface,'gray60',slider)
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
atoms = []
atoms.append(Component(39.9,0.34)) # Argon
atoms.append(Component(20.2,0.28)) # Neon
atoms.append(Component(83.8,0.40)) # Krypton
sliders = []
sliders.append(Slider('Argon Atoms','%4d','',N[0],0,100))  ## Change the minimum value of argon atoms to 0. Label the slider “Argon Atoms” ##
sliders.append(Slider('Neon Atoms','%4d','',N[1],0,100))   ## Added Neon Atoms and Slider ##
sliders.append(Slider('Krypton Atoms','%4d','',N[2],0,100))  ## Added Krypton Atoms and Slider ##
sliders.append(Slider('Temperature','%.1f','K',T,100,500))


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
init_pos = True
init_vel = True
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
                    if 'Atoms' in slider.lbl:
                        init_pos = True
                        init_vel = True
                    elif 'Temperature' in slider.lbl:
                        init_vel = True
        if event.type == pg.MOUSEBUTTONUP:
            alive = True
            for slider in sliders:
                slider.click = False

    for slider in sliders:
        slider.move((px,py))
        if 'Atoms' in slider.lbl:
            N[slider.idx] = round(slider.val)
        elif 'Temperature' in slider.lbl:
            T = slider.val

    ############################################################################
    # dynamics
    ############################################################################
    if alive:
        if init_pos:
            for iAtom,atom in enumerate(atoms):
                atom.init_pos(N[iAtom],box)
            init_pos = False
        if init_vel:
            for atom in atoms:
                atom.init_vel(T)
            init_vel = False
        for atom in atoms:
            atom.move(TIMESTEP,box)

    if N[0] + N[1] + N[2] < 2:
        for i in [0,1,2]:
            if sliders[i].val != 0:   ## Ensures atleast 1 atom is visually present, sometimes sliders get stuck in the 0 position ##
                sliders[i].click = True    ## or won't let the active atom slider move up, works best when sliders are slowly moved ##
                sliders[i].ddy = 0
                sliders[i].move((1,Slider.ymin))
                sliders[i].click = False

    Area = box[0]*box[1]  
    rho = sum(N[0:3])/Area
    Pressure = rho**(3/2)*kB*T  #units (g/mol/nm/ps2) ## Calculating Pressure using area number density and ideal gas law ##

    ############################################################################
    # render the screen
    ############################################################################
    surface.fill(SCREEN_COL)
    pg.draw.rect(surface,BOX_COL,(BOX_POS,BOX_DIM))
    for atom in atoms:
        atom.draw(surface)
    for slider in sliders:
        slider.draw(surface)
    font = pg.font.Font(None,3*12) # increased font size for better visibility #
    text = font.render('Pressure: %.4f bar' % (Pressure/k_to_u/mol_to_u*u_to_n*u_to_p**2/bar_to_Pa), True, 'white')  ## Internal Units to Bar ##
    textRect = text.get_rect()
    textRect.topleft = (BOX_POS[0]+10,BOX_POS[1]+10)
    surface.blit(text,textRect)
    pg.display.flip()

    clock.tick(FPS)

pg.display.quit()
pg.quit()
