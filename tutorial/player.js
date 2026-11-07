'use strict';
const video = document.querySelector('video');
const buttons = [...document.querySelectorAll('nav button')];
const speed = document.querySelector('#playback-speed');
speed.addEventListener('change', () => { video.playbackRate = Number(speed.value); });
video.addEventListener('ratechange', () => {
  const known = [...speed.options].find(option => Number(option.value) === video.playbackRate);
  if (known) speed.value = known.value;
});
function updateChapter() {
  let active = buttons[0];
  for (const button of buttons) {
    if (Number(button.dataset.start) <= video.currentTime) active = button;
  }
  for (const button of buttons) {
    button.classList.toggle('active', button === active);
    if (button === active) button.setAttribute('aria-current', 'true');
    else button.removeAttribute('aria-current');
  }
}
for (const button of buttons) {
  button.addEventListener('click', () => {
    video.currentTime = Number(button.dataset.start);
    video.play().catch(() => {});
    video.setAttribute('tabindex', '-1');
    video.focus({preventScroll: true});
    video.scrollIntoView({behavior: window.matchMedia('(prefers-reduced-motion: reduce)').matches ? 'auto' : 'smooth', block: 'start'});
    updateChapter();
  });
}
video.addEventListener('timeupdate', updateChapter);
video.addEventListener('loadedmetadata', updateChapter);
updateChapter();
