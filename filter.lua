function Image (img)
  img.src = img.src:gsub('^%.%./', './')
  return img
end
