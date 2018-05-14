# sudo apt-get install imagemagick

# convert -delay 20 -loop 0 *.jpg myimage.gif

# convert -resize 20% -delay 20 -loop 0 *.jpg myimage.gif

convert -resize 30% -delay 20 -loop 0 *.png test.gif

convert ns.png \( a.png le.png -append \) +append test.png

convert test.gif \( test_a.gif test_le.gif -append \) +append test_combine2.gif
