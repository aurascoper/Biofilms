/*
 * globe.js — the three.js scene, unchanged in behaviour from the pre-fusion page.
 *
 * Everything here is lifted from the original index.html with its comments
 * intact, because those comments record bugs that were already fixed once and
 * must not be reintroduced: the RA-as-hours error in raDecToVec3, the per-layer
 * camera-clamp swap, the north-to-south grid orientation.
 */

export const R_EARTH = 1.01;
export const R_CELESTIAL = 2.0;

export function latLonToVec3(lat, lon, r = 1.0) {
  const phi = (90 - lat) * Math.PI / 180;
  const theta = (lon + 180) * Math.PI / 180;
  return new THREE.Vector3(
    -r * Math.sin(phi) * Math.cos(theta),
    r * Math.cos(phi),
    r * Math.sin(phi) * Math.sin(theta));
}

// RA (degrees) + Dec (degrees) -> celestial shell.
// Delegating is the fix: the old body multiplied RA by 15, treating degrees as
// hours, and also flipped handedness relative to the terrestrial frame. Reusing
// latLonToVec3 corrects both at once and states the convention, which is that the
// celestial shell shares Earth's frame with RA read as longitude.
export function raDecToVec3(ra, dec, r = R_CELESTIAL) {
  return latLonToVec3(dec, ra, r);
}

export function createGlobe(container) {
  const scene = new THREE.Scene();
  const camera = new THREE.PerspectiveCamera(45, innerWidth / innerHeight, 0.1, 1000);
  camera.position.set(0, 0, 3.2);

  const renderer = new THREE.WebGLRenderer({ antialias: true });
  renderer.setSize(innerWidth, innerHeight);
  renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
  container.appendChild(renderer.domElement);

  const controls = new THREE.OrbitControls(camera, renderer.domElement);
  controls.enableDamping = true;
  controls.dampingFactor = 0.08;
  // One envelope. The old code swapped clamps when the stars layer was selected,
  // which is incoherent once Earth and sky are on screen together.
  controls.minDistance = 1.15;
  controls.maxDistance = 12;

  scene.add(new THREE.Points(
    new THREE.BufferGeometry().setFromPoints(
      Array.from({ length: 2000 }, () => new THREE.Vector3(
        (Math.random() - 0.5) * 200, (Math.random() - 0.5) * 200, (Math.random() - 0.5) * 200))),
    new THREE.PointsMaterial({ color: 0xffffff, size: 0.15, transparent: true, opacity: 0.6 })));

  // Earth is permanent now, never hidden.
  const earthGroup = new THREE.Group();
  // BARE_EARTH is the no-imagery colour. Held here because the imagery layer has
  // to restore it exactly when it detaches, and a texture on a near-black Phong
  // material renders as near-black.
  const earthMesh = new THREE.Mesh(
    new THREE.SphereGeometry(1, 64, 64),
    new THREE.MeshPhongMaterial({ color: 0x061219, transparent: true, opacity: 0.92 }));
  earthGroup.add(earthMesh);
  const lineMat = new THREE.LineBasicMaterial(
    { color: 0x11313d, transparent: true, opacity: 0.5 });
  for (let lat = -80; lat <= 80; lat += 20) {
    const pts = [];
    for (let lon = -180; lon <= 180; lon += 5) pts.push(latLonToVec3(lat, lon, 1.001));
    earthGroup.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(pts), lineMat));
  }
  for (let lon = -180; lon < 180; lon += 20) {
    const pts = [];
    for (let lat = -90; lat <= 90; lat += 5) pts.push(latLonToVec3(lat, lon, 1.001));
    earthGroup.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(pts), lineMat));
  }
  scene.add(earthGroup);

  const markerRoot = new THREE.Group();
  const bandRoot = new THREE.Group();
  scene.add(markerRoot, bandRoot);

  scene.add(new THREE.AmbientLight(0x404060, 0.6));
  const dir = new THREE.DirectionalLight(0xffffff, 1.0);
  dir.position.set(2, 1, 3);
  scene.add(dir);

  const raycaster = new THREE.Raycaster();
  raycaster.params.Line = { threshold: 0.01 };

  return { scene, camera, renderer, controls, markerRoot, bandRoot, raycaster,
           earthMesh, graticuleMaterial: lineMat, BARE_EARTH: 0x061219 };
}
