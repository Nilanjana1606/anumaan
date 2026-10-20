# Generates man/figures/logo.png. Run from the package root: Rscript data-raw/logo.R
library(ggplot2)
library(hexSticker)
sysfonts::font_add_google("Albert Sans", "albert", regular.wt = 600)
showtext::showtext_auto() # text size depends on dpi: p_size 12.5 is tuned for dpi = 300
out <- "man/figures/logo.png"

# Palette (Apple system colours): deep indigo body, light rim, orange accent
hex_fill <- "#2B2A6B"
hex_rim <- "#8E8CF0"
ink <- "#FFFFFF"
accent <- "#FF9F0A"

capsule <- function(w = 2.5, h = 1.05, n = 60) {
  r <- h / 2
  cx <- w / 2 - r
  a <- seq(-pi / 2, pi / 2, length.out = n)
  data.frame(
    x = c(cx + r * cos(a), -cx + r * cos(a + pi)),
    y = c(r * sin(a), r * sin(a + pi))
  )
}
rbar <- function(x0, h, w = 0.30, base = -0.32, n = 30) {
  r <- w / 2
  a <- seq(0, pi, length.out = n)
  data.frame(
    x = c(x0 + r, x0 + r * cos(a), x0 - r),
    y = c(base, base + h - r + r * sin(a), base)
  )
}
# upright bars whose bases follow the tilted body axis
tilt <- 22
ax <- function(u) c(u * cos(tilt * pi / 180), u * sin(tilt * pi / 180))
mk <- function(u, h, id) {
  p <- ax(u)
  cbind(rbar(p[1], h, base = p[2] - 0.27), id = id)
}
bars <- rbind(mk(-0.48, 0.30, 1), mk(0, 0.45, 2), mk(0.48, 0.60, 3))
tx <- seq(0, 1, length.out = 100)
tail <- data.frame(x = 1.25 + tx * 0.70, y = 0.10 * sin(tx * 3 * pi) * (1 - 0.3 * tx))
# slight tilt: reads as a microbe rather than a battery icon
rot <- function(d, deg = 22) {
  t <- deg * pi / 180
  transform(d, x = x * cos(t) - y * sin(t), y = x * sin(t) + y * cos(t))
}
tail <- rot(tail, tilt)
cap <- rot(capsule(), tilt)

g <- ggplot() +
  geom_path(data = tail, aes(x, y), colour = ink, linewidth = 1.4, lineend = "round") +
  geom_polygon(data = cap, aes(x, y), fill = NA, colour = ink, linewidth = 2.0) +
  geom_polygon(data = bars, aes(x, y, group = id, fill = id == 3), colour = NA) +
  scale_fill_manual(values = c(`FALSE` = ink, `TRUE` = accent), guide = "none") +
  coord_equal(xlim = c(-1.35, 1.95), ylim = c(-1.05, 1.15)) +
  theme_void() +
  theme_transparent()

sticker(g,
  package = "anumaan",
  p_family = "albert", p_fontface = "plain", p_size = 12.5, p_color = ink, p_y = 0.5,
  s_x = 1.04, s_y = 1.12, s_width = 1.45, s_height = 1.0,
  h_fill = hex_fill, h_color = hex_rim, h_size = 1.4,
  dpi = 300, filename = out
)

fav_dir <- "pkgdown/favicon"
dir.create(fav_dir, recursive = TRUE, showWarnings = FALSE)

pts <- function(d, k = 1) paste(sprintf("%.3f,%.3f", d$x, -d$y), collapse = " ")
glyph_svg <- paste0(
  '<g transform="translate(-0.30,0.05)">',
  sprintf('<polyline points="%s" fill="none" stroke="%s" stroke-width="0.17" stroke-linecap="round" stroke-linejoin="round"/>', pts(tail), ink),
  sprintf('<polygon points="%s" fill="none" stroke="%s" stroke-width="0.22" stroke-linejoin="round"/>', pts(cap), ink),
  paste(vapply(split(bars, bars$id), function(b) {
    sprintf(
      '<polygon points="%s" fill="%s"/>',
      pts(b), if (b$id[1] == 3) accent else ink
    )
  }, ""), collapse = ""),
  "</g>"
)

svg_doc <- function(body, view = "-2.2 -2.2 4.4 4.4") {
  sprintf('<svg xmlns="http://www.w3.org/2000/svg" viewBox="%s">%s</svg>', view, body)
}
a <- seq(90, 390, by = 60) * pi / 180 # pointy-top hexagon, as on the sticker
hex <- data.frame(x = 2.05 * cos(a), y = 2.05 * sin(a))
fav_svg <- svg_doc(paste0(
  sprintf(
    '<polygon points="%s" fill="%s" stroke="%s" stroke-width="0.16" stroke-linejoin="round"/>',
    pts(hex), hex_fill, hex_rim
  ),
  sprintf('<g transform="scale(1.05)">%s</g>', glyph_svg)
))
writeLines(fav_svg, file.path(fav_dir, "favicon.svg"))

touch_svg <- svg_doc(paste0(
  sprintf('<rect x="-2.2" y="-2.2" width="4.4" height="4.4" fill="%s"/>', hex_fill),
  sprintf('<g transform="scale(1.15)">%s</g>', glyph_svg)
))

render <- function(svg, px) magick::image_read_svg(svg, width = px, height = px)
magick::image_write(render(fav_svg, 96), file.path(fav_dir, "favicon-96x96.png"))
magick::image_write(render(touch_svg, 180), file.path(fav_dir, "apple-touch-icon.png"))
magick::image_write(magick::image_join(lapply(c(16, 32, 48), render, svg = fav_svg)),
  file.path(fav_dir, "favicon.ico"),
  format = "ico"
)
writeLines(sprintf('{
  "name": "anumaan",
  "short_name": "anumaan",
  "icons": [{ "src": "favicon-96x96.png", "sizes": "96x96", "type": "image/png" }],
  "theme_color": "%s",
  "background_color": "%s",
  "display": "standalone"
}', hex_fill, hex_fill), file.path(fav_dir, "site.webmanifest"))
