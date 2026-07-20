import pygame as pg
import numpy as np
rng = np.random.default_rng()

# box parms
bDim = pg.Vector2(500)

# screen parms
margin = 50
sDim = pg.Vector2(bDim.x,bDim.y+margin)
sColor = 'black'

# initialize pygame
pg.init()
surface = pg.display.set_mode(sDim)
clock = pg.time.Clock()
pg.display.set_caption('Pong')

# initial position of paddle
pPos = pg.Vector2(20,200)
pDim = pg.Vector2(10,100)
pVel = pg.Vector2(0,2)
pColor = 'blue'

# position and velocity of ball
R = 20
bPos = pg.Vector2(bDim.x-R,bDim.y/2)
bVel = pg.Vector2(2,rng.random()*2-1) 
bColor = 'red'

# scoreboard
frameColor = 'darkred'
scoreColor = 'white'
font = pg.font.Font(pg.font.get_default_font(),32)
hdr = font.render('Hits      High Score',True,scoreColor)
hdrRect = hdr.get_rect()
hdrRect.center = (250,bDim.y + margin//2)

# Power ups
powerups = []

running = True
alive = True
hit = False
score = 0

# load high score
def highscore_Load():
    try:
        with open('10_pong_highscore.txt','r') as f:
            hs = int(f.read())
    except:
        hs = 0
    return hs

highscore = highscore_Load()

extraPaddle = None

while running:

    # check to see if user shut window
    for event in pg.event.get():
        if event.type == pg.QUIT:
            running = False

    keys_pressed = pg.key.get_pressed()
    if keys_pressed[pg.K_q]:
        running = False

    if alive == True:
      pOld = pPos.copy() # store old position
      
      # paddle position
      if keys_pressed[pg.K_UP]:
          if pPos.y >= 0:
              pPos.y = max(0,pPos.y - pVel.y)
      if keys_pressed[pg.K_DOWN]:
          if pPos.y <= bDim.y:
              pPos.y = min(bDim.y - pDim.y,pPos.y + pVel.y)

      # ball position
      bPos += bVel

      # wall collisions
      if bPos.x >= bDim.x - R and bVel.x > 0:
          bVel.x *= -1
      if (bPos.y <= R and bVel.y < 0) or (bPos.y >= bDim.y - R and bVel.y > 0):
          bVel.y *= -1

      # paddle collision
      if bPos.x -R  <= pPos.x + pDim.x: 
          if bPos.y >= pPos.y and bPos.y <= pPos.y + pDim.y:
              hit = True
              score += 1
              bVel.x *= -1
              if pPos.y - pOld.y > 1:
                  bVel.y += pPos.y - pOld.y
              else:
                  bVel.y += rng.random()*2 -1 

              # Spawns powerup and randomly chooses type
              if score >= 3 and rng.random() < 0.25:
                  powerup_type = rng.choice(['green','white'])
                  y_pos = rng.uniform(0,bDim.y-20)
                  powerups.append({'pos':pg.Vector2(bDim.x-20, y_pos), 'type':powerup_type, 'color':'green' if powerup_type == 'green' else 'white'})
     
      # extra paddle collision
      if extraPaddle is not None:
          if bPos.x - R <= extraPaddle['pos'].x + extraPaddle['dim'].x:
              if bPos.y >= extraPaddle['pos'].y and bPos.y <= extraPaddle['pos'].y + extraPaddle['dim'].y:
                  hit = True
                  score += 1
                  bVel.x *= -1
                  if pPos.y - pOld.y > 1:
                      bVel.y += pPos.y - pOld.y
                  else:
                      bVel.y += rng.random()*2 -1 

                  extraPaddle['hits'] -= 1
                  if extraPaddle['hits'] <= 0:
                      extraPaddle = None    

      # powerups movement and collision  
      for powerup in powerups[:]:
            powerup['pos'].x -= 2  # move powerup left

            # power up collision with paddle
            if (powerup['pos'].x <= pPos.x + pDim.x and powerup['pos'].x + 20 >= pPos.x):
                    if (powerup['pos'].y + 10 >= pPos.y and powerup['pos'].y <= pPos.y + pDim.y):
                        if powerup['type'] == 'green':
                            bVel.x -= 1 # decreasing ball speed
                        
                        elif powerup['type'] == 'white':
                            rightmost_x = pPos.x + pDim.x
                            if extraPaddle is not None:
                                rightmost_x = max(rightmost_x, extraPaddle['pos'].x + extraPaddle['dim'].x)
                            extraPaddle = {'pos': pg.Vector2(rightmost_x + 50, pPos.y), 'dim': pg.Vector2(pDim.x, pDim.y),'hits':3}
                        powerups.remove(powerup)

            elif powerup['pos'].x + 20 < 0:
                powerups.remove(powerup)              


      # check for missed ball
      if bPos.x <= R:
            alive = False
            if score > highscore:
                highscore = score
                with open('10_pong_highscore.txt','w') as f:
                    f.write('%d' % highscore)

      # shorten paddle
      if hit and score%5 == 0:
          hit = False
          pDim.y = max(10,pDim.y - 10)

          if extraPaddle is not None:
              extraPaddle['dim'].y = pDim.y

          bVel.x += 1
          pVel.y += 1               

    # restart the game 
    if keys_pressed[pg.K_SPACE]:
        alive = True
        running = True
        hit = False
        score = 0
        pDim.y = 100
        pPos = pg.Vector2(20,200)
        bPos = pg.Vector2(bDim.x-R,bDim.y/2)
        bVel = pg.Vector2(2,rng.random()*2-1)
        powerups.clear()
        extraPaddle = None


    if alive == False:   
        GameOverTxt = font.render('Game Over', True, scoreColor)
        PressSpaceTxt = font.render('Press Space to Start', True, scoreColor)
        surface.blit(GameOverTxt, (50, 200))
        surface.blit(PressSpaceTxt, (50, 250))
        pg.display.flip()
        clock.tick(10) # text display keeps flickering, lowered fps to reduce flicker


    # update surface
    surface.fill(sColor)

    # render paddle
    pg.draw.rect(surface,pColor,(pPos,pDim))

    # render ball
    pg.draw.circle(surface,bColor,bPos,R)

    # render powerups 
    for powerup in powerups:
        pg.draw.rect(surface, powerup['color'], (powerup['pos'].x, powerup['pos'].y, 20, 20))

    if extraPaddle is not None:
        pg.draw.rect(surface, pColor, (extraPaddle['pos'], extraPaddle['dim']))

    # render score
    pg.draw.rect(surface,frameColor,(0,bDim.y,bDim.x,margin))
    surface.blit(hdr,hdrRect)
    scoreTxt = font.render('%d' % score, True, scoreColor)
    scoreRect = scoreTxt.get_rect()
    scoreRect.center = (200,bDim.y + margin//2)

    HighscoreTxt = font.render('%d' % highscore, True, scoreColor)
    HighscoreRect = HighscoreTxt.get_rect()
    HighscoreRect.center = (425,bDim.y + margin//2)
    surface.blit(HighscoreTxt,HighscoreRect)

    surface.blit(scoreTxt,scoreRect)
    pg.display.flip()

    # set fps
    clock.tick(60)
    
# exit out of pygame
pg.display.quit()
pg.quit()
