#!/usr/bin/env python
################################################################################
# IMPORT MODULES
################################################################################
import numpy as np
import pygame as pg
rng = np.random.default_rng()

################################################################################
# DEFINE CONSTANTS
################################################################################
# Game parameters
NAME = 'PacMan'
VERSION = 'Dev: Move Ghosts'

# margin parameters
BOTTOM_MARGIN, TOP_MARGIN = 0,0
LEFT_MARGIN, RIGHT_MARGIN = 0,0

# box parameters
BLOCK_SIZE = 100
BOX_POS = np.array([LEFT_MARGIN,TOP_MARGIN],dtype=int)
BOX_COL = 'black'

# screen parameters
SCREEN_COL = 'black'
FPS = 60

# key parameters
key_released = {}
key_released[pg.K_UP]    = True
key_released[pg.K_DOWN]  = True
key_released[pg.K_LEFT]  = True
key_released[pg.K_RIGHT] = True

# pacman parameters
RAD = BLOCK_SIZE // 2
VEL = BLOCK_SIZE // 10
DIR = {pg.K_UP:pg.Vector2(0,-1), pg.K_DOWN:pg.Vector2(0,1), pg.K_RIGHT:pg.Vector2(1,0), pg.K_LEFT:pg.Vector2(-1,0)}
COL = 'yellow'

# map parameters
MAP_FILE = 'map8.txt'
MAP_LIST = ['=',' ','p','.','g']
WALL_COL = 'darkblue'
PLT_COL = 'white'
PLT_SIZE = BLOCK_SIZE // 10

# ghost parameters
GST_COL = ['red','pink','cyan','orange']
GST_SIZE = BLOCK_SIZE // 2
GST_VEL = VEL
GST_MAX = 4


################################################################################
# DEFINE FUNCTIONS
################################################################################
def pbc(pos):
    margin = [LEFT_MARGIN, TOP_MARGIN]
    for d in range(2):
        if pos[d] > margin[d] + BOX_DIM[d]:
            pos[d] -= BOX_DIM[d]
        elif pos[d] < margin[d]:
            pos[d] += BOX_DIM[d]
    return pos

def read_map(inFilename, map_list):
    try:
        inFile = open(inFilename,'r')
    except IOError:
        nRows = 5
        nCols = 5
        world = [None]*nRows
        for i in range(nRows):
            world[i] = [' ']*nCols
        world[2][2] = 'p'
    else:
        nRows = 0
        nCols = 0
        world = []
        for line in inFile:
            row = []
            for block in line:
                if block in map_list:
                    row.append(block)
            iCols = len(row)
            if nCols == 0:
                nCols = iCols
            if nCols == iCols:
                world.append(row)
                nRows += 1
    shape = (nCols,nRows)
    # initialize objects
    walls = Walls(WALL_COL)
    pellets = Pellets(PLT_COL,PLT_SIZE)
    ghosts = []
    y = TOP_MARGIN
    for row in world:
        x = LEFT_MARGIN
        for block in row:
            if block == 'p':
                xc = x + BLOCK_SIZE//2
                yc = y + BLOCK_SIZE//2
                player = Pacman((xc,yc),COL)
            elif block == '=':
                walls.add_block((x,y))
            elif block == '.':
                xc = x + BLOCK_SIZE//2 - PLT_SIZE//2
                yc = y + BLOCK_SIZE//2 - PLT_SIZE//2
                pellets.add_pellet((xc,yc))
            elif block == 'g' and Ghosts.num < GST_MAX:
                ghosts.append(Ghosts((x,y),GST_SIZE))
            x += BLOCK_SIZE
        y += BLOCK_SIZE
            
    return shape,player,walls,pellets,ghosts

class Pacman(object):
    def __init__(self,pos,color):
        self.rad = BLOCK_SIZE//2
        self.pos = pos
        self.vel = pg.Vector2(0,0)
        self.new_vel = pg.Vector2(0,0)
        self.score = 0
        self.alive = True

    def move(self):
        # walls
        new = self.rect.move(self.new_vel)
        if new.collidelist(Walls.list) == -1:
            self.vel = self.new_vel
        old = self.rect.move(self.vel)
        if old.collidelist(Walls.list) == -1:
            self.pos = pbc(self.pos + self.vel)
        # pellets
        for pellet in Pellets.list:
            if self.rect.contains(pellet):
                self.score += 10
                Pellets.list.remove(pellet)
                Pellets.num -= 1
                if Pellets.num == 0:
                    self.alive = False
        # ghosts
        for ghost in ghosts:
            if self.rect.colliderect(ghost.small_rect):
                self.alive = False

    def draw(self,surface):
        self.rect = pg.draw.circle(surface,COL,self.pos,self.rad)

class Walls(object):
    
    num = 0
    list = []

    def __init__(self,color):
        self.color = color

    def add_block(self,pos):
        Walls.list.append(pg.Rect(pos,(BLOCK_SIZE,BLOCK_SIZE)))
        Walls.num += 1

    def draw(self,surface):
        for block in Walls.list:
            pg.draw.rect(surface,self.color,block)

class Pellets(object):

    num = 0
    list = []
    
    def __init__(self,color,size):
        self.color = color
        self.size = size
    
    def add_pellet(self,pos):
        Pellets.list.append(pg.Rect(pos,(self.size,self.size)))
        Pellets.num += 1

    def draw(self,surface):
        for pellet in Pellets.list:
            pg.draw.rect(surface,self.color,pellet)

class Ghosts(object):

    num = 0

    def __init__(self,pos,size):
        self.idx = Ghosts.num
        self.color = GST_COL[Ghosts.num%len(GST_COL)]
        self.size = size
        self.pos = pos
        self.big_rect = pg.Rect(self.pos,(BLOCK_SIZE,BLOCK_SIZE))
        self.vel = pg.Vector2(0,0)
        self.new_vel = GST_VEL * DIR[pg.K_LEFT]
        Ghosts.num += 1

    def move(self):
        new = self.big_rect.move(self.new_vel)
        if new.collidelist(Walls.list) == -1:
            self.vel = self.new_vel
        old = self.big_rect.move(self.vel)
        if old.collidelist(Walls.list) == -1:
            self.pos = pbc(self.pos + self.vel)
            self.big_rect = pg.Rect(self.pos,(BLOCK_SIZE,BLOCK_SIZE))
        else:
            self.new_vel = GST_VEL * DIR[rng.choice(list(DIR.keys()))]

    def draw(self,surface):
        x = self.pos[0] + BLOCK_SIZE//2 - self.size//2
        y = self.pos[1] + BLOCK_SIZE//2 - self.size//2
        self.small_rect = pg.Rect(x,y,self.size,self.size)
        pg.draw.rect(surface,self.color,self.small_rect)

################################################################################
# FILE MANAGEMENT
################################################################################
# check for required input files
shape,player,walls,pellets,ghosts = read_map(MAP_FILE,MAP_LIST)

# create output files

################################################################################
# MAIN CODE
################################################################################
# initialize box
BOX_DIM = np.array(shape,dtype=int)*BLOCK_SIZE
SCREEN_DIM = BOX_DIM + np.array([RIGHT_MARGIN+LEFT_MARGIN,BOTTOM_MARGIN+TOP_MARGIN])

# initial pygame
pg.init()
clock = pg.time.Clock()
surface = pg.display.set_mode(SCREEN_DIM)
pg.display.set_caption(f'{NAME:s} {VERSION:s}')
font12 = pg.font.Font(None,2*12)

# animation loop
running = True
while running:

    ############################################################################
    # input
    ############################################################################
    for event in pg.event.get():
        if event.type == pg.QUIT:
            running = False
        elif event.type == pg.KEYUP:
            key_released[event.key] = True

    key_pressed = pg.key.get_pressed()
    for key in DIR.keys():
        if key_pressed[key] and key_released[key] == True:
            key_released[key] = False
            player.new_vel = DIR[key]*VEL
        
    ############################################################################
    # render the screen
    ############################################################################
    surface.fill(SCREEN_COL)
    pg.draw.rect(surface,BOX_COL,(BOX_POS,BOX_DIM))
    # render map
    walls.draw(surface)
    pellets.draw(surface)
    # render pacman
    player.draw(surface)
    for ghost in ghosts:
        ghost.draw(surface)
    
    pg.display.flip()

    ############################################################################
    # dynamics
    ############################################################################
    if player.alive:
        player.move()
        for ghost in ghosts:
            ghost.move()
    
    if player.alive == False:
    #running == False:   
        PressSpaceTxt = font12.render('Press Space to Start', True, 'white')
        surface.blit(PressSpaceTxt, (500, 250))
        pg.display.flip()
        clock.tick(10) # text display keeps flickering, lowered fps to reduce flicker


        key_pressed = pg.key.get_pressed()
        if key_pressed[pg.K_SPACE]:
            alive = True
            pg.init()
            clock = pg.time.Clock()
            surface = pg.display.set_mode(SCREEN_DIM)
            pg.display.set_caption(f'{NAME:s} {VERSION:s}')

    clock.tick(FPS)

pg.display.quit()
pg.quit()
