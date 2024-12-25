p1_a = c(1.08, 1.99, 1.46, 1.21)
p1_b = c(1.48, 2.50, 2.62, 1.95)
p2_a = c(1.19, 2.10, 1.21, 0.96)
p2_b = c(0.62, 0.88, 0.68, 0.48)
p3_a = c(1.22, 1.91, 1.36, 0.90)
p3_b = c(0.65, 1.52, 1.32, 0.95)
p4_a = c(0.60, 1.10, 1.03, 0.61)
p4_b = c(0.32, 2.12, 1.48, 1.09)
p5_a = c(0.55, 1.00, 0.82, 0.52)
p5_b = c(1.48, 0.90, 0.75, 0.44)

t = c(1, 2, 3, 6)
plot(t, p1_a, type = "l", ylim = c(0, 3), ylab = "serum level", col = "purple")
lines(t, p1_b, col = "green")
lines(t, p2_a, col = "purple")
lines(t, p2_b, col = "green")
lines(t, p3_a, col = "purple")
lines(t, p3_b, col = "green")
lines(t, p4_a, col = "purple")
lines(t, p4_b, col = "green")
lines(t, p5_a, col = "purple")
lines(t, p5_b, col = "green")

# Begge medisinene virker som de øker antistoffet hos pasientene, men medisin A har en mer varierende effekt enn B virker det som.
# Det kan hende dette kommer av at denne medisinen er mer personavhengig.

# Vi har ikke hat noe om AUC og jeg har ikke boka =/