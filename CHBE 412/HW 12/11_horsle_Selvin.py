import numpy as np
import pygame as pg

# box parms
width = 2
bDim = pg.Vector2(50)
buffer = 10
nLetters = 5 
nAttempts = 6
boxColors = ['gray','yellow2','green4','black'] # wrong, close, right, default 

# screen parms
sDim = pg.Vector2(0)
sDim.x = nLetters*bDim.x + (nLetters+1)*buffer
sDim.y = (nAttempts+1)*bDim.y + (nAttempts+2)*buffer

# data structures
letters = -np.ones((nAttempts,nLetters),dtype=int)
status = -np.ones((nAttempts,nLetters),dtype=int)

#guesses
def get_solutions(filename):
    with open(filename,'r') as f:
        words = f.readlines()
    words = [w.strip() for w in words if len(w.strip()) == nLetters]
    return words

solutions_file = 'wordle_solutions.txt'
solutions = get_solutions(solutions_file)

# choose a word
#word = 'eying'
file = 'wordle_guesses.txt'

def get_random_word(filename):
    with open(file,'r') as f:
        words = f.readlines()
    words = [w.strip() for w in words if len(w.strip()) == nLetters]
    idx = np.random.randint(0,len(words))
    return words[idx]

word = get_random_word(file)

# pygame init
pg.init()
surface = pg.display.set_mode((sDim))
pg.display.set_caption('Horsle')
font = pg.font.Font(pg.font.get_default_font(),32)

# intialize counters & end conditions
Cele = ['Genius', 'Magnificent', 'Impressive', 'Splendid', 'Great', 'Phew']
message = ''
iLetters = 0
iAttempts = 0
nCorrect = 0
running = True
while running:
    # logic engine
    for event in pg.event.get():
        if event.type == pg.QUIT:
            running = False
        if event.type == pg.KEYUP:
            # alphabet
            if event.key >= pg.K_a and event.key <= pg.K_z and iLetters < 5:
                letters[iAttempts,iLetters] = event.key
                iLetters += 1
            # backspace
            if event.key == pg.K_BACKSPACE:
                iLetters = max(0,iLetters-1)
                letters[iAttempts,iLetters] = -1
            # return
            if event.key == pg.K_RETURN:
                nCorrect = 0

                if iLetters < 5:
                    message = 'Not enough letters'
                else:
                    guess = ''.join([pg.key.name(letters[iAttempts,i]) for i in range(nLetters)])

                    if guess in solutions:
                        for i in range(nLetters):
                            letter = pg.key.name(letters[iAttempts,i])
                            if letter == word[i]:
                                status[iAttempts,i] = 2
                                nCorrect += 1
                            elif letter in word:
                                status[iAttempts,i] = 1
                            else:
                                status[iAttempts,i] = 0
                    else:
                        message = 'Not in word list' ## Does not let me erase attempt afterwards ##
                    iAttempts += 1
                    iLetters = 0

        keys_pressed = pg.key.get_pressed()
        if iAttempts == nAttempts or nCorrect == nLetters:
            if keys_pressed[pg.K_SPACE]:
                word = get_random_word(file) # resets random word
                letters = -np.ones((nAttempts,nLetters),dtype=int)
                status = -np.ones((nAttempts,nLetters),dtype=int)
                iLetters = 0
                iAttempts = 0
                nCorrect = 0
            else: # ignores others inputs ## acts weird when person presses enter or backspace keys after game over ##
                continue
    
    # graphics
    surface.fill('black')
    Y = buffer
    for r in range(nAttempts):
        X = buffer
        for c in range(nLetters):
            pg.draw.rect(surface,'gray',(X,Y,bDim.x,bDim.y),width)
            pg.draw.rect(surface,boxColors[status[r,c]],(X+width,Y+width,bDim.x-2*width,bDim.y-2*width))
            if letters[r,c] != -1:
                letter = pg.key.name(letters[r,c])
                text = font.render(letter,True,'white')
                textRect = text.get_rect()
                textRect.center = (X + bDim.x//2,Y+bDim.y//2)
                surface.blit(text,textRect)
            X += bDim.x + buffer
        Y += bDim.y+buffer
    if nCorrect == nLetters:
        text = font.render(Cele[iAttempts-1],True,'white')
        textRect = text.get_rect()
        textRect.center = (sDim.x//2,Y+bDim.y//2)
        surface.blit(text,textRect)
    elif iAttempts == nAttempts:
        X = buffer
        for c in range(nLetters):
            text = font.render(word[c],True,'white')
            textRect = text.get_rect()
            textRect.center = (X+bDim.x//2,Y+bDim.y//2)
            surface.blit(text,textRect)
            X += bDim.x+buffer
    

    if iAttempts == nAttempts or nCorrect == nLetters:   
      #  surface.fill('black')
        pg.draw.rect(surface,'white', (15,90, 275, 100))
        font = pg.font.Font(pg.font.get_default_font(),20)
        GameOverTxt = font.render('Game Over', True, 'black')
        PressSpaceTxt = font.render('Press Space to Play Again', True, 'black')
        surface.blit(GameOverTxt, (25, 100))
        surface.blit(PressSpaceTxt, (25, 150))
        pg.display.flip()

    if message != '':
        font = pg.font.Font(pg.font.get_default_font(),20)
        text = font.render(message,True,'white')
        textRect = text.get_rect()
        textRect.center = (sDim.x//2,sDim.y - buffer - bDim.y//2)
        surface.blit(text,textRect)
        pg.display.flip()
        pg.time.delay(2000)
        message = ''


    pg.display.flip()
    
pg.display.quit()
pg.quit()
exit()
